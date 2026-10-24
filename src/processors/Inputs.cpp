#include "processors/Inputs.h"
#include "processors/WaveformProcessor.h"
#include "constructive_models/Magnetic.h"
#include "physical_models/MagnetizingInductance.h"
#include "support/Painter.h"

#include "support/Utils.h"
#include "support/NumberFormat.h"
#include "constructive_models/Coil.h"
#include "physical_models/WindingOhmicLosses.h"
#include "json.hpp"

#include <algorithm>
#include <cctype>
#include <cfloat>
#include <limits>
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
#include "support/Exceptions.h"


using json = nlohmann::json;

namespace OpenMagnetics {

// Cooley-Tukey FFT (in-place, breadth-first, decimation-in-frequency)
// Better optimized but less intuitive
// from https://rosettacode.org/wiki/Fast_Fourier_transform#C++
void fft(std::vector<std::complex<double>>& x) {
    // DFT
    unsigned int N = x.size(), k = N, n;
    double thetaT = std::numbers::pi / N;
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

// inverse fft (in-place)
void ifft(std::vector<std::complex<double>>& x)
{
    // conjugate the complex numbers
    for (size_t i = 0; i < x.size(); i++) {
        x[i] = std::conj(x[i]);
    }

    // forward fft
    fft( x );

    // conjugate the complex numbers again
    for (size_t i = 0; i < x.size(); i++) {
        x[i] = std::conj(x[i]);
    }

    // scale the numbers
    for (size_t i = 0; i < x.size(); i++) {
        x[i] /= x.size();
    }
}

double Inputs::calculate_waveform_average(Waveform waveform, std::optional<double> frequency) {
    // get_time() returns the optional BY VALUE, so indexing it inside the loop copied the whole
    // time axis twice per sample: O(N^2), minutes for an imported waveform with a few hundred
    // thousand points. Take it once.
    const auto time = waveform.get_time().value();
    const auto& data = waveform.get_data();
    if (frequency && WaveformProcessor::check_time_axis_spans_one_period(time, frequency.value())) {
        // Open axis (ABT #1713): sample i holds over [t_i, t_i+1), the last one over one more step
        // up to t0 + T. The trapezoid over t_0..t_N-1 left that last step out, divided by (N - 1) dt
        // and ramped every edge: an exactly balanced 688.5 V / -229.5 V rectangle (32 of 128
        // samples high) read a -1.807 V mean, which was then removed from the voltage before
        // integrating it and left a DC on its magnetizing current.
        const double lastStep = time[time.size() - 1] - time[time.size() - 2];
        double voltSeconds = data.back() * lastStep;
        for (size_t i = 0; i + 1 < data.size(); ++i) {
            voltSeconds += data[i] * (time[i + 1] - time[i]);
        }
        return voltSeconds / (time.back() + lastStep - time.front());
    }
    double integration = 0;
    double period = time.back() - time.front();
    for (size_t i = 0; i < data.size() - 1; ++i)
    {
        double area = (data[i + 1] + data[i]) / 2 * (time[i + 1] - time[i]);
        integration += area;
    }
    return integration / period;
}

Waveform Inputs::multiply_waveform(Waveform waveform, double scalarValue){
    auto& data = waveform.get_mutable_data();
    for (auto& datum : data) {
        datum *= scalarValue;
    }
    // A positive factor keeps the shape, so it keeps the label. Zero leaves an identically
    // zero waveform, which is no shape at all: a current scaled to zero kept "rectangular" on
    // 2 vertices and its duty was read from a third that does not exist (ABT #1878). A
    // negative factor flips the shape: a unipolar or flyback waveform is no longer one, and a
    // rectangle's duty becomes 1 - D, not the one its vertices are read for.
    if (!(scalarValue > 0)) {
        waveform.set_ancillary_label(std::nullopt);
    }
    return waveform;
}

Waveform Inputs::sum_waveform(Waveform waveform, double scalarValue){
    auto& data = waveform.get_mutable_data();
    for (auto& datum : data) {
        datum += scalarValue;
    }
    return waveform;
}

Waveform Inputs::add_waveforms(Waveform a, const Waveform& b) {
    auto& aData = a.get_mutable_data();
    const auto& bData = b.get_data();
    if (aData.size() != bData.size()) {
        throw std::invalid_argument("add_waveforms: size mismatch (" +
                                    std::to_string(aData.size()) + " vs " +
                                    std::to_string(bData.size()) + ")");
    }
    for (size_t i = 0; i < aData.size(); ++i) {
        aData[i] += bData[i];
    }
    return a;
}

Waveform Inputs::subtract_waveforms(Waveform a, const Waveform& b) {
    auto& aData = a.get_mutable_data();
    const auto& bData = b.get_data();
    if (aData.size() != bData.size()) {
        throw std::invalid_argument("subtract_waveforms: size mismatch (" +
                                    std::to_string(aData.size()) + " vs " +
                                    std::to_string(bData.size()) + ")");
    }
    for (size_t i = 0; i < aData.size(); ++i) {
        aData[i] -= bData[i];
    }
    return a;
}

Waveform Inputs::shift_waveform(Waveform waveform, int samples) {
    auto& data = waveform.get_mutable_data();
    const int n = static_cast<int>(data.size());
    if (n == 0) {
        throw std::invalid_argument("shift_waveform: empty waveform");
    }
    // Normalize shift into [0, n)
    int s = samples % n;
    if (s < 0) s += n;
    if (s == 0) return waveform;
    std::vector<double> shifted(n);
    for (int i = 0; i < n; ++i) {
        shifted[(i + s) % n] = data[i];
    }
    data = shifted;
    return waveform;
}

SignalDescriptor Inputs::add_signals(const SignalDescriptor& a, const SignalDescriptor& b, double frequency) {
    if (!a.get_waveform() || !b.get_waveform()) {
        throw std::invalid_argument("add_signals: both signals must have waveforms");
    }
    auto wa = a.get_waveform().value();
    auto wb = b.get_waveform().value();
    auto sa = calculate_sampled_waveform(wa, frequency);
    auto sb = calculate_sampled_waveform(wb, frequency);
    auto sum = add_waveforms(sa, sb);
    SignalDescriptor out;
    out.set_waveform(sum);
    out.set_harmonics(calculate_harmonics_data(sum, frequency));
    out.set_processed(calculate_processed_data(out, sum, true));
    return out;
}

SignalDescriptor Inputs::subtract_signals(const SignalDescriptor& a, const SignalDescriptor& b, double frequency) {
    if (!a.get_waveform() || !b.get_waveform()) {
        throw std::invalid_argument("subtract_signals: both signals must have waveforms");
    }
    auto wa = a.get_waveform().value();
    auto wb = b.get_waveform().value();
    auto sa = calculate_sampled_waveform(wa, frequency);
    auto sb = calculate_sampled_waveform(wb, frequency);
    auto diff = subtract_waveforms(sa, sb);
    SignalDescriptor out;
    out.set_waveform(diff);
    out.set_harmonics(calculate_harmonics_data(diff, frequency));
    out.set_processed(calculate_processed_data(out, diff, true));
    return out;
}

SignalDescriptor Inputs::shift_signal(const SignalDescriptor& s, int samples, double frequency) {
    if (!s.get_waveform()) {
        throw std::invalid_argument("shift_signal: signal must have a waveform");
    }
    auto sampled = calculate_sampled_waveform(s.get_waveform().value(), frequency);
    auto shifted = shift_waveform(sampled, samples);
    SignalDescriptor out;
    out.set_waveform(shifted);
    out.set_harmonics(calculate_harmonics_data(shifted, frequency));
    out.set_processed(calculate_processed_data(out, shifted, true));
    return out;
}

bool is_close_enough(double x, double y, double error){
    return fabs(x - y) <= error;
}

// Every winding's current on ONE time grid, and at which instants it carries current (ABT #1467).
// The old classifier compared the windings' RAW waveforms index by index whenever a winding's
// current and voltage had the same number of points, ignoring their time vectors: a 3-winding
// flyback whose secondaries were tabulated on another grid looked simultaneously conducting and
// lost its magnetizing DC. Here every current is resampled onto the same numberPoints over the
// period first. Conduction is judged on the current (Ampere: a winding carries ampere-turns iff it
// carries current), against 1 % of that winding's own peak, the classifier's existing figure.
// std::nullopt when no winding has a current: the operating point is given by its voltages
// alone (transformer mode), and the currents say nothing about the magnetizing MMF. When some
// windings have a current and another has none, the ampere-turns cannot be summed: that throws
// (ABT #1814), where it used to return std::nullopt and the magnetizing current lost its DC
// without a word.
std::optional<Inputs::WindingConduction> Inputs::calculate_winding_conduction(const OperatingPoint& operatingPoint, std::optional<size_t> numberPoints) {
    const auto excitations = operatingPoint.get_excitations_per_winding();
    if (excitations.empty()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "calculate_winding_conduction: the operating point has no excitations");
    }
    std::optional<size_t> windingWithCurrent;
    std::optional<size_t> windingWithoutCurrent;
    for (size_t windingIndex = 0; windingIndex < excitations.size(); ++windingIndex) {
        if (excitations[windingIndex].get_current()) {
            if (!windingWithCurrent) {
                windingWithCurrent = windingIndex;
            }
        }
        else if (!windingWithoutCurrent) {
            windingWithoutCurrent = windingIndex;
        }
    }
    if (!windingWithCurrent) {
        return std::nullopt;
    }
    if (windingWithoutCurrent) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "calculate_winding_conduction: winding " + std::to_string(windingWithoutCurrent.value()) +
                                    " has no current while winding " + std::to_string(windingWithCurrent.value()) +
                                    " has one: the magnetizing MMF is the sum of every winding's ampere-turns, so the DC of the magnetizing current"
                                    " cannot be told from the currents unless every winding carries one");
    }
    size_t gridPoints = settings.get_inputs_number_points_sampled_waveforms();
    for (size_t windingIndex = 0; windingIndex < excitations.size(); ++windingIndex) {
        const auto current = excitations[windingIndex].get_current();
        if (!current->get_waveform()) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "calculate_winding_conduction: winding " + std::to_string(windingIndex) +
                                        " has a current without a waveform; processed values carry no timing to compare the windings on");
        }
        size_t windingPoints = 1;
        while (windingPoints < current->get_waveform()->get_data().size()) {
            windingPoints *= 2;
        }
        gridPoints = std::max(gridPoints, windingPoints);
        if (excitations[windingIndex].get_frequency() != excitations[0].get_frequency()) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "calculate_winding_conduction: winding " + std::to_string(windingIndex) + " is excited at " +
                                        std::to_string(excitations[windingIndex].get_frequency()) + " Hz, winding 0 at " + std::to_string(excitations[0].get_frequency()) + " Hz");
        }
    }
    if (numberPoints) {
        if (numberPoints.value() == 0 || (numberPoints.value() & (numberPoints.value() - 1)) != 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "calculate_winding_conduction: the grid must have a power-of-2 number of points, got " + std::to_string(numberPoints.value()));
        }
        gridPoints = numberPoints.value();
    }

    WindingConduction conduction;
    for (size_t windingIndex = 0; windingIndex < excitations.size(); ++windingIndex) {
        auto sampled = calculate_sampled_waveform(excitations[windingIndex].get_current()->get_waveform().value(),
                                                  excitations[windingIndex].get_frequency(), gridPoints, gridPoints);
        auto data = sampled.get_data();
        if (data.size() != gridPoints) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "calculate_winding_conduction: winding " + std::to_string(windingIndex) + " sampled to " +
                                        std::to_string(data.size()) + " points instead of " + std::to_string(gridPoints));
        }
        double peak = 0;
        for (auto value : data) {
            peak = std::max(peak, std::fabs(value));
        }
        std::vector<bool> conducting(gridPoints);
        for (size_t pointIndex = 0; pointIndex < gridPoints; ++pointIndex) {
            conducting[pointIndex] = std::fabs(data[pointIndex]) > peak * 0.01;
        }
        conduction.currents.push_back(data);
        conduction.conducting.push_back(conducting);
    }
    return conduction;
}

std::vector<bool> Inputs::calculate_magnetizing_conduction_instants(const WindingConduction& conduction, size_t windingIndex) {
    // The instants at which winding windingIndex's current is the magnetizing current (referred to
    // it): it conducts, and the primary (winding 0) hands it the magnetizing MMF rather than
    // carrying load current with it. For the primary: no other winding conducts. For any other
    // winding: the primary does not. A flyback's primary on-time and its secondaries' off-time
    // qualify; so does a forward's reset interval (the primary alone, its current returning the
    // magnetizing energy), while a forward's on-time, where primary and secondary carry the load
    // current together, does not.
    const size_t numberWindings = conduction.conducting.size();
    if (windingIndex >= numberWindings) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "calculate_magnetizing_conduction_instants: winding " + std::to_string(windingIndex) +
                                    " of an operating point with " + std::to_string(numberWindings) + " windings");
    }
    const size_t gridPoints = conduction.conducting[windingIndex].size();
    std::vector<bool> instants(gridPoints, false);
    for (size_t pointIndex = 0; pointIndex < gridPoints; ++pointIndex) {
        if (!conduction.conducting[windingIndex][pointIndex]) {
            continue;
        }
        bool alone = true;
        if (windingIndex == 0) {
            for (size_t otherIndex = 1; otherIndex < numberWindings; ++otherIndex) {
                if (conduction.conducting[otherIndex][pointIndex]) {
                    alone = false;
                    break;
                }
            }
        }
        else {
            alone = !conduction.conducting[0][pointIndex];
        }
        instants[pointIndex] = alone;
    }
    return instants;
}

bool windings_alternate(const Inputs::WindingConduction& conduction) {
    // Windings that do not always conduct together: an adjacent pair differs on more than 10 % of
    // the period (the classifier's existing figure).
    const auto& conducting = conduction.conducting;
    const size_t gridPoints = conducting[0].size();
    size_t numberPointsNotConductingAtTheSameTime = 0;
    for (size_t pointIndex = 0; pointIndex < gridPoints; ++pointIndex) {
        for (size_t windingIndex = 0; windingIndex + 1 < conducting.size(); ++windingIndex) {
            if (conducting[windingIndex][pointIndex] != conducting[windingIndex + 1][pointIndex]) {
                numberPointsNotConductingAtTheSameTime++;
                break;
            }
        }
    }
    return numberPointsNotConductingAtTheSameTime > gridPoints * 0.1;
}

bool Inputs::include_dc_offset_into_magnetizing_current(OperatingPoint operatingPoint, std::vector<double> turnsRatios) {
    // Determine if this is truly a single-winding component (inductor) by checking both
    // the turns ratios and the actual number of excitations in the operating point
    bool onlyOneWinding = turnsRatios.size() == 0 && operatingPoint.get_excitations_per_winding().size() <= 1;
    if (onlyOneWinding) {
        return true;
    }

    // Windings that do not always conduct together (a flyback, or a forward's reset interval):
    // somewhere one of them carries the magnetizing MMF alone, so the currents give the DC
    // (calculate_magnetizing_current_dc_anchor). Windings that always conduct together (a
    // transformer) carry load current on top of it, and their DC is not anchored. Compared on one
    // time grid (ABT #1467); a pair differing on more than 10 % of the period counts, as before.
    // The voltage's zero intervals used to veto the DC too (include_dc_offset_into_magnetizing_current_rosano):
    // a DCM flyback's dead time is 0 V, so its magnetizing current lost its DC once it came from
    // Faraday instead of the FLYBACK label's triangle (ABT #1467).
    auto conduction = calculate_winding_conduction(operatingPoint);
    if (!conduction) {
        return false;
    }
    return windings_alternate(conduction.value());
}

double Inputs::calculate_magnetizing_current_dc_anchor(const OperatingPoint& operatingPoint, size_t windingIndex, Waveform voltageSampledWaveform, double magnetizingInductance) {
    // The DC of winding windingIndex's magnetizing current, from the currents (ABT #907, #1467).
    // Faraday gives the magnetizing current's shape, i_F(t) = (1/L) int v dt with zero mean; the DC
    // is where the currents say what the magnetizing MMF is: at the instants where this winding
    // carries it alone (calculate_magnetizing_conduction_instants), by Ampere i_m = i_winding.
    //  - No winding conducts somewhere (a dead time): i_m = 0 there. Preferred, as it does not
    //    depend on the L_m the currents were simulated with.
    //  - The winding's current peaks there (a flyback): peak(i_m) = peak(i_winding), anchored on
    //    the current's exact peak and the excursion's maximum (ABT #907: the commutation steps of
    //    the current move no flux, so its midpoint or pk-pk do not describe i_m).
    //  - It peaks while sharing the load current (a forward's on-time): the peak is load current
    //    and anchors nothing. The largest current of the magnetizing-only instants (the start of
    //    the forward's reset) anchors instead, at the same instant of i_F.
    //  - It never carries the magnetizing MMF alone (a forward's secondary): no DC.
    const auto excitations = operatingPoint.get_excitations_per_winding();
    if (windingIndex >= excitations.size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "calculate_magnetizing_current_dc_anchor: winding " + std::to_string(windingIndex) +
                                    " of an operating point with " + std::to_string(excitations.size()) + " windings");
    }
    if (magnetizingInductance <= 0) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "magnetizingInductance cannot be zero or negative");
    }
    auto excitation = excitations[windingIndex];
    if (!excitation.get_current()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "calculate_magnetizing_current_dc_anchor: winding " + std::to_string(windingIndex) +
                                    " has no current to anchor the DC of its magnetizing current on");
    }
    // Every winding has a current here, or calculate_winding_conduction throws (ABT #1814).
    auto conduction = calculate_winding_conduction(operatingPoint).value();
    return calculate_magnetizing_current_dc_anchor(excitation, conduction, windingIndex, voltageSampledWaveform, magnetizingInductance);
}

double Inputs::calculate_magnetizing_current_dc_anchor(OperatingPointExcitation excitation, const WindingConduction& conduction, size_t windingIndex,
                                                       Waveform voltageSampledWaveform, double magnetizingInductance) {
    if (!excitation.get_current()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "calculate_magnetizing_current_dc_anchor: winding " + std::to_string(windingIndex) +
                                    " has no current to anchor the DC of its magnetizing current on");
    }
    if (windingIndex >= conduction.currents.size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "calculate_magnetizing_current_dc_anchor: winding " + std::to_string(windingIndex) +
                                    " of a conduction map with " + std::to_string(conduction.currents.size()) + " windings");
    }
    const auto& current = conduction.currents[windingIndex];

    // Instants where no winding conducts and the winding voltage is zero (a forward's dead time
    // after its reset, a DCM flyback's dead time): every winding's ampere-turns are zero, so by
    // Ampere the magnetizing current is zero there, whatever inductance the currents were
    // simulated with. The voltage must be zero too: a winding voltage means the flux is still
    // moving, so currents that read zero there are a gap in the data, not a dead time (the
    // mas-flyback-3-winding-lv fixture's secondaries start conducting 2.8 us after the primary
    // stops, under the full -8 V reset voltage). This anchor needs only the
    // voltage's shape, so it is preferred to the ones below, which convert a current into flux
    // through this core's L_m and inherit any mismatch between that L_m and the one behind the
    // currents (the forward-228w fixture's reset current, 0.55 A, is a ~1.2 mH reset; on a core with
    // another L_m, anchoring on it shifted B by 0.14 T).
    auto voltageOnGrid = calculate_sampled_waveform(voltageSampledWaveform, excitation.get_frequency(), current.size(), current.size()).get_data();
    if (voltageOnGrid.size() != current.size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "calculate_magnetizing_current_dc_anchor: the voltage sampled to " +
                                    std::to_string(voltageOnGrid.size()) + " points, the currents to " + std::to_string(current.size()));
    }
    double voltagePeak = 0;
    for (auto value : voltageOnGrid) {
        voltagePeak = std::max(voltagePeak, std::fabs(value));
    }
    std::vector<size_t> deadInstants;
    for (size_t pointIndex = 0; pointIndex < current.size(); ++pointIndex) {
        if (std::fabs(voltageOnGrid[pointIndex]) > voltagePeak * 0.01) {
            continue;
        }
        bool anyConducting = false;
        for (const auto& windingConducting : conduction.conducting) {
            if (windingConducting[pointIndex]) {
                anyConducting = true;
                break;
            }
        }
        if (!anyConducting) {
            deadInstants.push_back(pointIndex);
        }
    }
    if (!deadInstants.empty()) {
        // The same mean-free integral the magnetizing current is built from (calculate_magnetizing_current).
        double meanVoltage = calculate_steady_state_mean_voltage(excitation, voltageSampledWaveform, "Excitation '" + excitation.get_name().value_or("unnamed") + "'");
        auto centeredMagnetizingWaveform = calculate_integral_waveform(sum_waveform(voltageSampledWaveform, -meanVoltage), true);
        centeredMagnetizingWaveform = multiply_waveform(centeredMagnetizingWaveform, 1.0 / magnetizingInductance);
        auto excursionOnGrid = calculate_sampled_waveform(centeredMagnetizingWaveform, excitation.get_frequency(), current.size(), current.size()).get_data();
        if (excursionOnGrid.size() != current.size()) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "calculate_magnetizing_current_dc_anchor: the volt-second integral sampled to " +
                                        std::to_string(excursionOnGrid.size()) + " points, the currents to " + std::to_string(current.size()));
        }
        double excursionWhileDead = 0;
        for (auto pointIndex : deadInstants) {
            excursionWhileDead += excursionOnGrid[pointIndex];
        }
        return -excursionWhileDead / deadInstants.size();
    }

    auto instants = calculate_magnetizing_conduction_instants(conduction, windingIndex);
    std::optional<size_t> aloneMaximumIndex;
    std::optional<double> sharedMaximum;
    for (size_t pointIndex = 0; pointIndex < current.size(); ++pointIndex) {
        if (instants[pointIndex]) {
            if (!aloneMaximumIndex || current[pointIndex] > current[aloneMaximumIndex.value()]) {
                aloneMaximumIndex = pointIndex;
            }
        }
        else if (conduction.conducting[windingIndex][pointIndex]) {
            if (!sharedMaximum || current[pointIndex] > sharedMaximum.value()) {
                sharedMaximum = current[pointIndex];
            }
        }
    }
    if (!aloneMaximumIndex) {
        return 0;
    }

    // The same mean-free integral the magnetizing current is built from (calculate_magnetizing_current).
    double meanVoltage = calculate_steady_state_mean_voltage(excitation, voltageSampledWaveform, "Excitation '" + excitation.get_name().value_or("unnamed") + "'");
    auto centeredMagnetizingWaveform = calculate_integral_waveform(sum_waveform(voltageSampledWaveform, -meanVoltage), true);
    centeredMagnetizingWaveform = multiply_waveform(centeredMagnetizingWaveform, 1.0 / magnetizingInductance);

    if (!sharedMaximum || current[aloneMaximumIndex.value()] >= sharedMaximum.value()) {
        auto currentSignal = excitation.get_current().value();
        if (!currentSignal.get_processed() || !currentSignal.get_processed()->get_peak()) {
            auto currentWaveform = currentSignal.get_waveform().value();
            currentSignal.set_processed(calculate_processed_data(currentWaveform, excitation.get_frequency()));
        }
        const auto& excursion = centeredMagnetizingWaveform.get_data();
        double excursionAboveMean = *max_element(excursion.begin(), excursion.end());
        return currentSignal.get_processed()->get_peak().value() - excursionAboveMean;
    }
    auto excursionOnGrid = calculate_sampled_waveform(centeredMagnetizingWaveform, excitation.get_frequency(), current.size(), current.size()).get_data();
    if (excursionOnGrid.size() != current.size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "calculate_magnetizing_current_dc_anchor: the volt-second integral sampled to " +
                                    std::to_string(excursionOnGrid.size()) + " points, the currents to " + std::to_string(current.size()));
    }
    return current[aloneMaximumIndex.value()] - excursionOnGrid[aloneMaximumIndex.value()];
}


double Inputs::try_guess_duty_cycle(Waveform waveform, WaveformLabel label, double frequency) {
    // Delegates to WaveformProcessor, which fixed the same classifier bug
    // (ABT #602): three unconditional `return 0.5` defaults for genuinely
    // unrecognised shapes, and a second-difference heuristic that divided an
    // index into `data` by numberPointsSampledWaveforms (a 128-vs-512
    // mismatch for imported waveforms) and thresholded at 5% of the maximum
    // (duty = 1 for any DC-biased signal). This used to be a verbatim,
    // independently-bugged twin of that logic — kept as a delegating wrapper
    // so the two can no longer drift apart, matching every other Inputs::X
    // waveform-DSP method in this file.
    return WaveformProcessor::try_guess_duty_cycle(waveform, label, frequency,
                                                   settings.get_inputs_number_points_sampled_waveforms());
}

bool Inputs::is_standardized(SignalDescriptor signal) {
    if (!signal.get_waveform()) {
        return false;
    }

    if (signal.get_waveform() && !signal.get_waveform()->get_time()) {
        return false;
    }

    return true;
}

// In case the waveform comes defined with processed data only, we create a MAS format waveform from it, as the rest of
// the code depends on it.
SignalDescriptor Inputs::standardize_waveform(SignalDescriptor signal, double frequency) {
    // SignalDescriptor standardized_signal = signal;
    SignalDescriptor standardized_signal(signal);
    if (!signal.get_waveform()) {
        if (!signal.get_processed() && !signal.get_harmonics()) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Signal is not processed");
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

Waveform Inputs::reconstruct_signal(Harmonics harmonics, double frequency) {
    size_t numberPoints = std::max(settings.get_inputs_number_points_sampled_waveforms(), 16 * round_up_size_to_power_of_2(harmonics.get_frequencies().back() / frequency));
    std::vector<double> data = std::vector<double>(numberPoints, 0);


    for (size_t harmonicIndex = 0; harmonicIndex < harmonics.get_frequencies().size(); ++harmonicIndex) {
        auto amplitude = harmonics.get_amplitudes()[harmonicIndex];
        auto frequencyMultiplier = harmonics.get_frequencies()[harmonicIndex] / frequency;
        // Use numberPoints (not numberPoints-1) to match FFT convention:
        // N samples span [0, T) with step T/N, period boundary is NOT included.
        auto totalAngle = 2 * std::numbers::pi / numberPoints * frequencyMultiplier;
        for (size_t i = 0; i < numberPoints; ++i) {
            if (harmonics.get_frequencies()[harmonicIndex] > 0) {
                double angle = i * totalAngle;
                data[i] += sin(angle) * amplitude;
            }
            else {
                data[i] += amplitude;
            }
        }
    }
    // Generate time array matching FFT convention: N samples spanning [0, T)
    // with step T/N, where T = 1/frequency. The period endpoint is NOT included.
    double period = 1. / roundFloat(frequency, 9);
    std::vector<double> time(numberPoints);
    for (size_t i = 0; i < numberPoints; ++i) {
        time[i] = i * period / numberPoints;
    }
    Waveform waveform;
    waveform.set_data(data);
    waveform.set_time(time);
    return waveform;
}

void Inputs::throw_if_processed_cannot_define_waveform(const ProcessedWaveform& processed, const std::string& what) {
    auto label = processed.get_label();
    std::string labelName(magic_enum::enum_name(label));
    if (label == WaveformLabel::CUSTOM || label == WaveformLabel::RECTANGULAR_DCM) {
        throw InvalidInputException(ErrorCode::MISSING_DATA,
            "The " + what + " has no waveform and its processed label " + labelName +
            " has no parametric form to build one from: supply the sampled waveform");
    }
    std::vector<std::string> missing;
    if (!processed.get_peak_to_peak()) {
        missing.push_back("peakToPeak");
    }
    if (label != WaveformLabel::SINUSOIDAL && !processed.get_duty_cycle()) {
        missing.push_back("dutyCycle");
    }
    const bool withDeadTime = label == WaveformLabel::TRIANGULAR_WITH_DEADTIME ||
                              label == WaveformLabel::RECTANGULAR_WITH_DEADTIME ||
                              label == WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME ||
                              label == WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME;
    if (withDeadTime && !processed.get_dead_time()) {
        missing.push_back("deadTime");
    }
    if (!missing.empty()) {
        std::string fields;
        for (auto& field : missing) {
            fields += (fields.empty() ? "" : ", ") + field;
        }
        throw InvalidInputException(ErrorCode::MISSING_DATA,
            "The " + what + " has no waveform, and its processed data cannot define a " + labelName +
            " one: missing " + fields);
    }
}

Waveform Inputs::create_waveform(ProcessedWaveform processed, double frequency) {
    return WaveformProcessor::create_waveform(processed, frequency,
                                              settings.get_inputs_number_points_sampled_waveforms());
}


Waveform Inputs::create_waveform(WaveformLabel label, double peakToPeak, double frequency, double dutyCycle, double offset, double deadTime, double skew, double phase) {
    return WaveformProcessor::create_waveform(label, peakToPeak, frequency, dutyCycle, offset, deadTime, skew, phase,
                                              settings.get_inputs_number_points_sampled_waveforms());
}

bool Inputs::is_waveform_sampled(Waveform waveform) {
    if (!waveform.get_time()) {
        return false;
    }
    else if (is_waveform_imported(waveform)) {
        return is_size_power_of_2(waveform.get_data());
    }
    else {
        return waveform.get_data().size() == settings.get_inputs_number_points_sampled_waveforms();
    }

}

bool Inputs::is_waveform_imported(Waveform waveform) {
    if (!waveform.get_time()) {
        return false;
    }
    else {
        return waveform.get_data().size() > settings.get_inputs_number_points_sampled_waveforms();
    }

}

bool Inputs::has_flyback_label(const OperatingPointExcitation& excitation) {
    auto current = excitation.get_current();
    if (!current) {
        return false;
    }
    if (current->get_waveform() && current->get_waveform()->get_ancillary_label()) {
        WaveformLabel ancillaryLabel = current->get_waveform()->get_ancillary_label().value();
        if (ancillaryLabel == WaveformLabel::FLYBACK_PRIMARY || ancillaryLabel == WaveformLabel::FLYBACK_SECONDARY) {
            return true;
        }
    }
    if (current->get_processed()) {
        auto label = current->get_processed()->get_label();
        if (label == WaveformLabel::FLYBACK_PRIMARY || label == WaveformLabel::FLYBACK_SECONDARY) {
            return true;
        }
    }
    return false;
}

bool Inputs::is_multiport_inductor(OperatingPoint operatingPoint, std::optional<std::vector<IsolationSide>> isolationSides) {
    auto excitations = operatingPoint.get_excitations_per_winding();
    if (excitations.size() < 2) {
        return false;
    }
    else {
        if (isolationSides) {
            bool allSameIsolationSide = true;
            auto isolationSidesValue = isolationSides.value();
            for (size_t windingIndex = 1; windingIndex < isolationSidesValue.size(); ++windingIndex) {
                if (isolationSidesValue[0] != isolationSidesValue[windingIndex]) {
                    allSameIsolationSide = false;
                }
            }
            if (allSameIsolationSide) {
                return true;
            }
        }

        OperatingPointExcitation excitation = Inputs::get_primary_excitation(operatingPoint);
        bool hasFlybackLabel = has_flyback_label(excitation);

        if (hasFlybackLabel) {
            // Guard: If another winding shares the primary's isolation side, this is a
            // forward converter with a demagnetization winding, not a flyback/coupled inductor.
            // Forward converters should use the voltage integration path instead.
            if (isolationSides && isolationSides->size() >= 2) {
                auto primarySide = isolationSides.value()[0];
                for (size_t windingIndex = 1; windingIndex < isolationSides->size(); ++windingIndex) {
                    if (isolationSides.value()[windingIndex] == primarySide) {
                        return false;
                    }
                }
            }
            return true;
        }
    }
    return false;
}


bool Inputs::can_be_common_mode_choke(OperatingPoint operatingPoint, std::optional<std::vector<IsolationSide>> isolationSides) {
    auto excitations = operatingPoint.get_excitations_per_winding();
    // CMC topologies: 2 (line + neutral), 3 (three-phase delta), or 4
    // (three-phase + neutral). Four-winding was previously excluded, which
    // forced those designs through the multiport-inductor path and pumped
    // the full line current into the core-flux calculation.
    if (excitations.size() < 2 || excitations.size() > 4) {
        return false;
    }
    // MKF signs a winding's ampere-turns by its isolation side: windings on different sides oppose,
    // so the same current on a primary and a secondary cancels in the core. That is a transformer
    // carrying load current, not common-mode current; reading it as a CMC stored the load current
    // as the magnetizing current (ABT #1467: a 1:1 transformer's 3 A instead of V / (w L_m) = 0.16 A).
    if (isolationSides) {
        if (isolationSides->size() < excitations.size()) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "can_be_common_mode_choke: " + std::to_string(isolationSides->size()) +
                                        " isolation sides for " + std::to_string(excitations.size()) + " windings");
        }
        for (size_t windingIndex = 1; windingIndex < excitations.size(); ++windingIndex) {
            if (isolationSides.value()[windingIndex] != isolationSides->front()) {
                return false;
            }
        }
    }

    if (!operatingPoint.get_excitations_per_winding()[0].get_current()) {
        throw std::invalid_argument("Current is missing");
    }
    auto primaryCurrent = operatingPoint.get_excitations_per_winding()[0].get_current().value();

    if (primaryCurrent.get_harmonics()) {
        // Harmonics path: all winding harmonics must match (same CM source seen by every winding).
        // Analytical CMC produces exact matches; simulated CMC may have tiny numerical discrepancies
        // from the ngspice circuit (parasitic-cap asymmetry). If all harmonics match exactly, return
        // true immediately. If any harmonic mismatches, fall through to the processed-data RMS check
        // below — which handles simulation-derived operating points where harmonics differ slightly
        // but the underlying physical current is the same CM noise signal.
        bool harmonicsMatch = true;
        for (size_t windingIndexIndex = 1; windingIndexIndex < excitations.size() && harmonicsMatch; ++windingIndexIndex) {
            auto secondaryCurrent = excitations[windingIndexIndex].get_current().value();

            if (primaryCurrent.get_harmonics() && secondaryCurrent.get_harmonics()) {
                auto primaryHarmonics = primaryCurrent.get_harmonics().value();
                auto secondaryHarmonics = secondaryCurrent.get_harmonics().value();
                for (size_t harmonicIndex = 0; harmonicIndex < primaryHarmonics.get_frequencies().size(); ++harmonicIndex) {
                    if (!is_close_enough(primaryHarmonics.get_frequencies()[harmonicIndex], secondaryHarmonics.get_frequencies()[harmonicIndex], 0.0001)) {
                        harmonicsMatch = false;
                        break;
                    }
                    if (!is_close_enough(primaryHarmonics.get_amplitudes()[harmonicIndex], secondaryHarmonics.get_amplitudes()[harmonicIndex], 0.0001)) {
                        harmonicsMatch = false;
                        break;
                    }
                }
            }
        }
        if (harmonicsMatch) return true;
        // Fall through to processed-data check for simulation-derived operating points.
    }

    // Time-domain / simulation path: simulate_realistic_cmc sets both harmonics and processed data.
    // The harmonics may not match exactly due to numerical simulation noise (parasitic capacitance
    // asymmetry in the ngspice model). Check instead that all winding RMS amplitudes agree within
    // 5%, which is sufficient to distinguish a symmetric CM choke from an asymmetric inductor.
    if (primaryCurrent.get_processed() && primaryCurrent.get_processed()->get_rms()) {
        double primaryRms = primaryCurrent.get_processed()->get_rms().value();
        for (size_t windingIndex = 1; windingIndex < excitations.size(); ++windingIndex) {
            auto secondaryCurrent = excitations[windingIndex].get_current();
            if (!secondaryCurrent || !secondaryCurrent->get_processed() || !secondaryCurrent->get_processed()->get_rms()) {
                return false;
            }
            double secondaryRms = secondaryCurrent->get_processed()->get_rms().value();
            // is_close_enough takes an absolute tolerance; the documented intent is 5% relative
            if (!is_close_enough(primaryRms, secondaryRms, 0.05 * primaryRms)) {
                return false;
            }
        }
        return true;
    }

    return false;
}

bool Inputs::is_multiport_inductor_magnetizing_current_from_voltage(const OperatingPointExcitation& primaryExcitation) {
    return has_flyback_label(primaryExcitation) && primaryExcitation.get_voltage().has_value();
}

SignalDescriptor Inputs::get_multiport_inductor_magnetizing_current(OperatingPoint operatingPoint,
                                                                   double magnetizingInductance,
                                                                   const std::vector<double>& turnsRatios,
                                                                   std::optional<std::vector<IsolationSide>> isolationSides) {
    // ABT #1467: this used to rebuild a TRIANGULAR current from the primary current's offset, peak
    // to peak and duty cycle, rising over D and falling over the whole of 1 - D. A DCM flyback's
    // dead time was lost (its flux, and so its core losses, read from the wrong shape), and for
    // windings on one side the other windings' ampere-turns were ignored.
    OperatingPointExcitation excitation = Inputs::get_primary_excitation(operatingPoint);
    if (!excitation.get_current()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Multiport inductor: the primary excitation has no current");
    }
    SignalDescriptor magnetizingCurrent;
    size_t numberWindings = operatingPoint.get_excitations_per_winding().size();
    if (is_multiport_inductor_magnetizing_current_from_voltage(excitation)) {
        // The windings conduct alternately: the flux comes from the primary voltage (Faraday),
        // anchored where the primary carries the whole magnetizing MMF alone.
        magnetizingCurrent = calculate_magnetizing_current(operatingPoint, 0, magnetizingInductance, turnsRatios, false);
    }
    else {
        // No primary voltage on a flyback, or every winding on one isolation side: the ampere-turn
        // sum of the winding currents, referred to the primary, each winding signed by its isolation
        // side as calculate_ampere_turn_current does. By Ampere that sum IS the magnetizing MMF only
        // while the given currents are the windings' true currents: for a flyback it is exact while
        // either winding conducts and through the dead time (both zero), but it carries no flux the
        // currents do not show (no ringing, no leakage-commutation overlap), which Faraday would.
        if (!isolationSides || isolationSides->empty()) {
            throw InvalidInputException(ErrorCode::MISSING_DATA,
                "Multiport inductor: without a flyback primary voltage, the magnetizing current is the windings' ampere-turn sum, which needs their isolation sides");
        }
        std::vector<IsolationSide> sides;
        for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
            sides.push_back(windingIndex < isolationSides->size() ? isolationSides.value()[windingIndex] : get_isolation_side_from_index(windingIndex));
        }
        if (!has_flyback_label(excitation)) {
            for (auto side : sides) {
                if (side != sides.front()) {
                    throw InvalidInputException(ErrorCode::INVALID_INPUT,
                        "Multiport inductor without a flyback label: its windings must share one isolation side for their ampere-turns to add");
                }
            }
            // calculate_ampere_turn_current needs a primary-side reference; with every winding on
            // one side, which side it is does not change the sum.
            std::fill(sides.begin(), sides.end(), IsolationSide::PRIMARY);
        }
        auto ampereTurnCurrent = calculate_ampere_turn_current(operatingPoint, turnsRatios, sides);
        auto waveform = ampereTurnCurrent.get_waveform().value();
        magnetizingCurrent.set_waveform(waveform);
        magnetizingCurrent.set_harmonics(calculate_harmonics_data(waveform, waveform, excitation.get_frequency()));
    }
    auto sampledWaveform = magnetizingCurrent.get_waveform().value();
    if (!is_size_power_of_2(sampledWaveform.get_data())) {
        sampledWaveform = Inputs::calculate_sampled_waveform(sampledWaveform, excitation.get_frequency());
        magnetizingCurrent.set_waveform(sampledWaveform);
    }
    magnetizingCurrent.set_processed(calculate_processed_data(magnetizingCurrent, sampledWaveform, true));
    return magnetizingCurrent;
}

SignalDescriptor Inputs::get_common_mode_choke_magnetizing_current(OperatingPoint operatingPoint) {
    OperatingPointExcitation excitation = Inputs::get_primary_excitation(operatingPoint);

    // ── Physics recap for CMCs ───────────────────────────────────────────
    // Each winding carries: DM line current + CM noise current.
    // DM currents in the windings are opposite in real life (line vs neutral)
    // but MKF stores every winding's current "into the dot", so after that
    // sign flip the DM shows up as the difference and the CM as the sum.
    // By CMC symmetry all windings carry identical waveforms, so:
    //   sum  = 2·I_cm   (drives core flux)
    //   diff = 0        (DM, cancels — not a core excitation)
    double frequency = Inputs::get_switching_frequency(excitation);

    // The core's MMF is the sum of the windings' ampere-turns. A CMC's windings have the same turns and their
    // currents are stored into the dot, so the current that drives the core is the windings' mean: their
    // common-mode part is the same in all of them and stays, their differential-mode part cancels. Taking the
    // primary's own swing read a differential-mode line point (+/-14 A on the two windings, opposite, which
    // the harmonic amplitudes cannot tell from identical currents) as 14 A of common-mode current: 18.8 T in
    // a 4-turn R10KZ toroid whose real flux there is ~0 (ABT #2022). For a symmetric CMC point the mean's
    // swing is the primary's, as before.
    // We model the magnetizing current as a triangular of the mean's p-p centered around zero (no DC flux
    // bias in a CMC).
    std::vector<double> meanCurrent;
    std::optional<Waveform> meanWaveform;
    const auto& excitations = operatingPoint.get_excitations_per_winding();
    for (const auto& windingExcitation : excitations) {
        if (!windingExcitation.get_current()) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Common-mode choke: a winding carries no current, so the windings' mean current cannot be formed");
        }
        // Each current is sampled over the period it is given over, its excitation's (a line-cycle record
        // spans the line period, not the switching period the triangle below is built at).
        auto current = windingExcitation.get_current().value();
        if (!current.get_waveform()) {
            current = Inputs::standardize_waveform(current, windingExcitation.get_frequency());
        }
        auto sampled = Inputs::calculate_sampled_waveform(current.get_waveform().value(), windingExcitation.get_frequency());
        if (!meanWaveform) {
            meanWaveform = sampled;
            meanCurrent.assign(sampled.get_data().size(), 0.0);
        }
        else if (sampled.get_data().size() != meanCurrent.size()) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Common-mode choke: the windings' currents sample to different lengths (" +
                                        std::to_string(meanCurrent.size()) + " and " + std::to_string(sampled.get_data().size()) + ")");
        }
        for (size_t index = 0; index < meanCurrent.size(); ++index) {
            meanCurrent[index] += sampled.get_data()[index] / excitations.size();
        }
    }
    double peakToPeak = *std::max_element(meanCurrent.begin(), meanCurrent.end()) - *std::min_element(meanCurrent.begin(), meanCurrent.end());

    ProcessedWaveform triangularProcessed;
    triangularProcessed.set_label(WaveformLabel::TRIANGULAR);
    // Centered around zero, as documented above: a CMC carries no DC flux bias, so the
    // triangular swings ±pp/2. The previous offset of pp/2 shifted the waveform to 0..pp,
    // doubling the peak the saturation/loss models saw — contradicting this very comment.
    triangularProcessed.set_offset(0);
    triangularProcessed.set_peak_to_peak(peakToPeak);
    auto waveform = create_waveform(triangularProcessed, frequency);
    SignalDescriptor magnetizingCurrent;
    auto sampledWaveform = Inputs::calculate_sampled_waveform(waveform, frequency);
    magnetizingCurrent.set_waveform(sampledWaveform);
    magnetizingCurrent.set_harmonics(calculate_harmonics_data(waveform, sampledWaveform, frequency));
    magnetizingCurrent.set_processed(calculate_processed_data(magnetizingCurrent, sampledWaveform, true));

    return magnetizingCurrent;
}

bool Inputs::can_be_differential_mode_choke(OperatingPoint operatingPoint, bool isDmcTopology) {
    // DMC dispatch fires when:
    //   1) caller has identified the topology as DMC
    //      (designRequirements.topology == DIFFERENTIAL_MODE_CHOKE), AND
    //   2) the operating point has 2, 3, or 4 windings — implying
    //      SINGLE_PHASE_BALANCED, THREE_PHASE, or THREE_PHASE_WITH_NEUTRAL.
    //
    // Single-phase 1-winding DMC (1 winding) is intentionally NOT
    // intercepted here — it's a plain series inductor and the existing
    // numberWindings==1 path in MagnetizingInductance.cpp handles it
    // correctly (full line current with DC bias).
    if (!isDmcTopology) {
        return false;
    }
    auto excitations = operatingPoint.get_excitations_per_winding();
    return excitations.size() >= 2 && excitations.size() <= 4;
}

SignalDescriptor Inputs::get_differential_mode_choke_magnetizing_current(OperatingPoint operatingPoint) {
    // ── DMC effective core current synthesis ─────────────────────────────
    // MKF convention: every winding's current is stored "into the dot",
    // matching DifferentialModeChoke::process_operating_points which uses
    // the same phase sign on every winding (the physical reversal of
    // winding sense is already absorbed into the sign convention). The
    // flux-driving "magnetizing" current the core sees is therefore the
    // point-by-point SUM of all winding currents for every DMC config:
    //
    //   2 windings (SINGLE_PHASE_BALANCED): i_L + i_N (both same sign
    //     into-the-dot → 2·I_DM driving flux additively).
    //
    //   3 windings (THREE_PHASE 3-wire):    i_a + i_b + i_c (vector sum;
    //     in a perfectly balanced load this is zero — no excitation —
    //     and imbalance / zero-sequence is what we filter).
    //
    //   4 windings (THREE_PHASE_WITH_NEUTRAL): i_a + i_b + i_c + i_n.
    //
    // 1-winding DMC (SINGLE_PHASE) is handled by the regular inductor
    // path and never lands here (can_be_differential_mode_choke returns
    // false for 1 winding).
    OperatingPointExcitation primaryExcitation = Inputs::get_primary_excitation(operatingPoint);
    double frequency = primaryExcitation.get_frequency();

    auto excitations = operatingPoint.get_excitations_per_winding();
    if (excitations.size() < 2 || excitations.size() > 4) {
        throw std::invalid_argument("DMC magnetizing current synthesis requires 2, 3, or 4 windings, got " +
                                    std::to_string(excitations.size()));
    }

    std::vector<Waveform> windingCurrents;
    windingCurrents.reserve(excitations.size());
    for (size_t i = 0; i < excitations.size(); ++i) {
        if (!excitations[i].get_current()) {
            throw std::invalid_argument("DMC: winding " + std::to_string(i) + " has no current");
        }
        auto current = excitations[i].get_current().value();
        if (!current.get_waveform()) {
            throw std::invalid_argument("DMC: winding " + std::to_string(i) + " current has no waveform");
        }
        auto wf = current.get_waveform().value();
        if (wf.get_data().empty() || !is_size_power_of_2(wf.get_data())) {
            wf = calculate_sampled_waveform(wf, frequency);
        }
        windingCurrents.push_back(wf);
    }

    size_t N = windingCurrents[0].get_data().size();
    for (size_t i = 1; i < windingCurrents.size(); ++i) {
        if (windingCurrents[i].get_data().size() != N) {
            throw std::invalid_argument("DMC: per-winding waveform sizes mismatch — winding 0 has " +
                std::to_string(N) + ", winding " + std::to_string(i) + " has " +
                std::to_string(windingCurrents[i].get_data().size()));
        }
    }

    std::vector<double> combined(N, 0.0);
    for (size_t s = 0; s < N; ++s) {
        double sum = 0.0;
        for (size_t w = 0; w < windingCurrents.size(); ++w) {
            sum += windingCurrents[w].get_data()[s];
        }
        combined[s] = sum;
    }

    Waveform combinedWaveform;
    combinedWaveform.set_data(combined);
    if (windingCurrents[0].get_time()) {
        combinedWaveform.set_time(windingCurrents[0].get_time().value());
    }

    SignalDescriptor magnetizingCurrent;
    magnetizingCurrent.set_waveform(combinedWaveform);
    magnetizingCurrent.set_harmonics(calculate_harmonics_data(combinedWaveform, frequency));
    magnetizingCurrent.set_processed(calculate_processed_data(magnetizingCurrent, combinedWaveform, true));

    return magnetizingCurrent;
}

SignalDescriptor Inputs::calculate_ampere_turn_current(const OperatingPoint& operatingPoint,
                                                       const std::vector<double>& turnsRatios,
                                                       const std::vector<IsolationSide>& isolationSides) {
    const auto& excitations = operatingPoint.get_excitations_per_winding();
    const size_t numberWindings = excitations.size();
    if (numberWindings == 0) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "calculate_ampere_turn_current: the operating point has no excitations");
    }
    if (isolationSides.size() != numberWindings) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "calculate_ampere_turn_current: " + std::to_string(isolationSides.size()) + " isolation sides for " + std::to_string(numberWindings) + " windings; the sign of each winding's contribution depends on its side");
    }
    if (!turnsRatios.empty() && turnsRatios.size() != numberWindings - 1) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "calculate_ampere_turn_current: " + std::to_string(turnsRatios.size()) + " turns ratios for " + std::to_string(numberWindings) + " windings (expected one per winding after the first)");
    }
    // N_k / N_0: winding 0 is 1, winding k is 1 / (Np/Nk).
    std::vector<double> turnsAgainstFirst(numberWindings, 1.0);
    for (size_t k = 1; k < numberWindings && !turnsRatios.empty(); ++k) {
        if (turnsRatios[k - 1] <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "calculate_ampere_turn_current: turns ratio " + std::to_string(turnsRatios[k - 1]) + " is not positive");
        }
        turnsAgainstFirst[k] = 1.0 / turnsRatios[k - 1];
    }
    size_t reference = numberWindings;
    for (size_t k = 0; k < numberWindings; ++k) {
        if (isolationSides[k] == IsolationSide::PRIMARY) {
            reference = k;
            break;
        }
    }
    if (reference == numberWindings) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "calculate_ampere_turn_current: no winding is on the primary side, so there is no reference winding");
    }

    std::optional<double> frequency;
    std::vector<double> sum;
    std::vector<double> time;
    for (size_t k = 0; k < numberWindings; ++k) {
        const auto excitation = excitations[k];
        const auto current = excitation.get_current();
        if (!current || !current->get_waveform()) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "calculate_ampere_turn_current: winding " + std::to_string(k) + " has no current waveform; processed values carry no phase to align the windings with");
        }
        if (frequency && std::fabs(excitation.get_frequency() - frequency.value()) > 1e-9 * frequency.value()) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "calculate_ampere_turn_current: winding " + std::to_string(k) + " is excited at " + std::to_string(excitation.get_frequency()) + " Hz, winding 0 at " + std::to_string(frequency.value()) + " Hz");
        }
        frequency = excitation.get_frequency();
        auto sampled = calculate_sampled_waveform(current->get_waveform().value(), frequency.value());
        const auto& data = sampled.get_data();
        const double scale = (isolationSides[k] == IsolationSide::PRIMARY ? 1.0 : -1.0) * turnsAgainstFirst[k] / turnsAgainstFirst[reference];
        if (sum.empty()) {
            sum.assign(data.size(), 0.0);
            time = sampled.get_time().value();
        }
        if (data.size() != sum.size()) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "calculate_ampere_turn_current: winding " + std::to_string(k) + " sampled to " + std::to_string(data.size()) + " points, winding 0 to " + std::to_string(sum.size()));
        }
        for (size_t i = 0; i < data.size(); ++i) {
            sum[i] += scale * data[i];
        }
    }
    Waveform waveform;
    waveform.set_data(sum);
    waveform.set_time(time);
    SignalDescriptor ampereTurnCurrent;
    ampereTurnCurrent.set_waveform(waveform);
    ampereTurnCurrent.set_processed(calculate_processed_data(waveform, frequency.value()));
    return ampereTurnCurrent;
}

Waveform Inputs::calculate_sampled_waveform(Waveform waveform, double frequency, std::optional<size_t> numberPoints, std::optional<size_t> maximumNumberPoints) {
    return WaveformProcessor::calculate_sampled_waveform(waveform, frequency, numberPoints,
                                                         settings.get_inputs_number_points_sampled_waveforms(),
                                                         maximumNumberPoints);
}

double Inputs::calculate_max_volt_seconds(const OperatingPointExcitation& excitation) {
    if (!excitation.get_voltage()) {
        return 0.0;
    }
    // MAS::OperatingPointExcitation::get_voltage() (and every other MAS-
    // generated optional getter on this path) returns
    // `std::optional<...>` BY VALUE (MAS.hpp:4650). Binding
    //     const auto& voltage = excitation.get_voltage().value();
    // takes a reference to the contents of the temporary optional, which
    // dies at the end of the full-expression — leaving `voltage` dangling
    // and silently returning an empty optional from
    // `voltage.get_waveform()` even when the original excitation has a
    // populated voltage waveform. Symptom upstream was the transformer
    // saturation-aware initial-N (CoreAdviserDataset.cpp:552) computing
    // `nFromSat = 0`, falling back to 5 turns, and then the saturation
    // filter rejecting every ferrite core for PSFB / SRC / LLC / Vienna
    // designs because B_peak = V·s / (N·Ae) ≫ B_sat at N=5.
    // Copy each level by value so all references stay valid for the
    // whole function.
    const auto voltage = excitation.get_voltage().value();
    const auto waveformOpt = voltage.get_waveform();
    if (waveformOpt.has_value() && waveformOpt.value().get_time().has_value()) {
        const auto& wf = waveformOpt.value();
        const auto& data = wf.get_data();
        const auto timeOpt = wf.get_time();
        const auto& time = timeOpt.value();
        const size_t n = std::min(data.size(), time.size());
        if (n < 2) return 0.0;

        // Limit integration to at most ONE switching period. Multi-cycle
        // waveforms (typically ngspice transient output) with a small
        // numerical DC bias on V accumulate spurious drift in the running
        // integral: max|∫V dt| grows linearly with the simulation duration
        // instead of reflecting the per-cycle flux excursion. The relevant
        // quantity for saturation sizing is the per-period peak flux
        // ΔΨ = max|∫₀^T V dt| (Faraday), which is well-defined on one
        // period regardless of how many periods the caller supplied.
        //
        // If frequency is unknown (legacy path), fall back to integrating
        // the whole waveform — the historical behaviour.
        const double frequency = excitation.get_frequency();
        const double tStart = time.front();
        const double tEnd = (frequency > 0.0)
                                ? tStart + 1.0 / frequency
                                : time.back() + 1.0;  // never reached
        double integral = 0.0;
        double maxVoltSeconds = 0.0;
        for (size_t j = 0; j + 1 < n; ++j) {
            if (time[j] >= tEnd) break;
            integral += data[j] * (time[j + 1] - time[j]);
            maxVoltSeconds = std::max(maxVoltSeconds, std::abs(integral));
        }
        return maxVoltSeconds;
    }
    // Sinusoidal fallback from processed peak: V·t_max = V_peak / ω.
    const auto processedOpt = voltage.get_processed();
    if (processedOpt.has_value() && processedOpt.value().get_peak().has_value()) {
        // DOCUMENTED DEFAULT: with no usable frequency this rough estimator
        // assumes 100 kHz (typical SMPS); the result is only a sizing heuristic
        // and callers needing exactness must provide the frequency.
        double frequency = excitation.get_frequency() > 0 ? excitation.get_frequency() : 100000.0;
        double omega = 2.0 * std::numbers::pi * frequency;
        if (omega > 0) {
            return processedOpt.value().get_peak().value() / omega;
        }
    }
    return 0.0;
}

double Inputs::calculate_max_volt_seconds(const OperatingPoint& operatingPoint) {
    double maxVoltSeconds = 0.0;
    for (const auto& exc : operatingPoint.get_excitations_per_winding()) {
        maxVoltSeconds = std::max(maxVoltSeconds, calculate_max_volt_seconds(exc));
    }
    return maxVoltSeconds;
}

SignalDescriptor Inputs::calculate_induced_voltage(OperatingPointExcitation& excitation, double magnetizingInductance, bool compress) {
    if (!excitation.get_current()->get_waveform()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Current waveform is missing");
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
        // FIXED: BUG-04 — Use period-based wrapping
        double wrapPeriod = time.back() - time.front();
        source.push_back(source[0]);
        time.push_back(time[time.size() - 1] + wrapPeriod);
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
    auto sampledWaveform = Inputs::calculate_sampled_waveform(voltageWaveform, excitation.get_frequency());
    voltageSignalDescriptor.set_harmonics(calculate_harmonics_data(voltageWaveform, sampledWaveform, excitation.get_frequency()));
    voltageSignalDescriptor.set_processed(calculate_processed_data(voltageSignalDescriptor, sampledWaveform, true));
    if (!compress) {
        voltageSignalDescriptor.set_waveform(sampledWaveform);
    }

    resultWaveform.set_data(derivative);
    return voltageSignalDescriptor;
}

Waveform Inputs::calculate_derivative_waveform(Waveform waveform) {
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
        // FIXED: BUG-04 — Use period-based wrapping
        double wrapPeriod = tempTime.back() - tempTime.front();
        tempData.push_back(tempData[0]);
        tempTime.push_back(tempTime[tempTime.size() - 1] + wrapPeriod);
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

Waveform Inputs::calculate_integral_waveform(Waveform waveform, bool subtractAverage) {
    std::vector<double> data = waveform.get_data();
    std::vector<double> time = waveform.get_time().value();
    std::vector<double> integration;
    Waveform resultWaveform(waveform);
    // The integral is a different shape from its integrand: a rectangular voltage integrates
    // to a triangular magnetizing current. Copying the voltage's label onto it made a 3- or
    // 4-vertex triangle "rectangular" or "bipolarRectangular", whose offset and duty were then
    // read at vertex indices it does not have (ABT #1878). The integral's shape is the one its
    // data shows, so it carries no label.
    resultWaveform.set_ancillary_label(std::nullopt);

    // FIXED: BUG-02 — Ensure integration result matches time vector size
    double integral = 0;
    size_t integrationSize = std::min(data.size(), time.size());
    integration.push_back(integral);
    for (size_t i = 0; i < integrationSize - 1; ++i) {
        double timePerPoint = time[i + 1] - time[i];
        integral += data[i] * timePerPoint;
        integration.push_back(integral);
    }
    resultWaveform.set_data(integration);

    if (subtractAverage) {
        // The integration constant is the integral's average over the whole PERIOD (ABT #1713).
        // Sample i holds data[i] over [t_i, t_i+1); the last sample holds over one more step,
        // [t_N-1, t_N-1 + dt), which closes the period of a sampled waveform (N points spanning
        // [0, T)). The trapezoid of the integral over t_0..t_N-1 alone left that interval out and
        // divided by (N - 1) dt: a zero-volt-second voltage kept ~1 % of its peak as DC.
        if (integrationSize < 2) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "calculate_integral_waveform: the waveform needs at least 2 samples to remove the integration constant");
        }
        double lastStep = time[integrationSize - 1] - time[integrationSize - 2];
        if (!(lastStep > 0)) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "calculate_integral_waveform: the last two samples share an instant, so the hold of the last sample (and the period) is undefined");
        }
        double closedIntegral = integration.back() + data[integrationSize - 1] * lastStep;
        double area = 0;
        for (size_t i = 0; i + 1 < integration.size(); ++i) {
            area += (integration[i] + integration[i + 1]) / 2 * (time[i + 1] - time[i]);
        }
        area += (integration.back() + closedIntegral) / 2 * lastStep;
        double period = time[integrationSize - 1] + lastStep - time[0];
        resultWaveform = sum_waveform(resultWaveform, -area / period);
    }

    // Taken once: get_time() copies the whole axis on every call (see calculate_waveform_average).
    const auto resultTime = resultWaveform.get_time().value();
    const auto& resultData = resultWaveform.get_data();
    std::vector<double> distinctData;
    std::vector<double> distinctTime;
    for (size_t i = 0; i < resultData.size(); ++i){
        if (distinctData.size() != 0) {
            if (resultData[i] == distinctData.back() && resultTime[i] == distinctTime.back())
                continue;
        }

        distinctData.push_back(resultData[i]);
        distinctTime.push_back(resultTime[i]);
    }

    resultWaveform.set_data(distinctData);
    resultWaveform.set_time(distinctTime);

    return resultWaveform;
}

SignalDescriptor Inputs::add_offset_to_excitation(SignalDescriptor signalDescriptor,
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
    auto sampledWaveform = Inputs::calculate_sampled_waveform(waveform, frequency);
    signalDescriptor.set_harmonics(calculate_harmonics_data(waveform, sampledWaveform, frequency));
    signalDescriptor.set_processed(calculate_processed_data(signalDescriptor, sampledWaveform, true, signalDescriptor.get_processed()));
    return signalDescriptor;
}

OperatingPointExcitation Inputs::get_excitation_with_proportional_current(OperatingPointExcitation excitation, double proportion) {
    if (!excitation.get_current()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Excitation is missing current");
    }
    if (!excitation.get_current()->get_waveform()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Current is missing waveform");
    }
    auto current = excitation.get_current().value();
    auto multipliedWaveform = multiply_waveform(current.get_waveform().value(), proportion);
    current.set_waveform(multipliedWaveform);
    auto sampledCurrentWaveform = calculate_sampled_waveform(multipliedWaveform, excitation.get_frequency());
    current.set_harmonics(calculate_harmonics_data(multipliedWaveform, sampledCurrentWaveform, excitation.get_frequency()));
    // The basic processed data (peak, peak to peak, offset) scale with the waveform, so they are calculated again:
    // handing over the previous processed data would keep its values. Only the label survives, as the shape does,
    // which is only for a positive proportion (see multiply_waveform).
    auto previousProcessedCurrent = current.get_processed();
    auto processedCurrent = calculate_processed_data(current, sampledCurrentWaveform, true, std::nullopt);
    if (previousProcessedCurrent && proportion > 0) {
        processedCurrent.set_label(previousProcessedCurrent->get_label());
    }
    current.set_processed(processedCurrent);

    excitation.set_current(current);
    return excitation;
}

OperatingPointExcitation Inputs::get_excitation_with_proportional_voltage(OperatingPointExcitation excitation, double proportion) {
    if (!excitation.get_voltage()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Excitation is missing voltage");
    }
    if (!excitation.get_voltage()->get_waveform()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "voltage is missing waveform");
    }
    auto voltage = excitation.get_voltage().value();
    auto multipliedWaveform = multiply_waveform(voltage.get_waveform().value(), proportion);
    voltage.set_waveform(multipliedWaveform);
    auto sampledVoltageWaveform = calculate_sampled_waveform(multipliedWaveform, excitation.get_frequency());
    voltage.set_harmonics(calculate_harmonics_data(multipliedWaveform, sampledVoltageWaveform, excitation.get_frequency()));
    // The basic processed data (peak, peak to peak, offset) scale with the waveform, so they are calculated again:
    // handing over the previous processed data would keep its values. Only the label survives, as the shape does,
    // which is only for a positive proportion (see multiply_waveform).
    auto previousProcessedVoltage = voltage.get_processed();
    auto processedVoltage = calculate_processed_data(voltage, sampledVoltageWaveform, true, std::nullopt);
    if (previousProcessedVoltage && proportion > 0) {
        processedVoltage.set_label(previousProcessedVoltage->get_label());
    }
    voltage.set_processed(processedVoltage);
    excitation.set_voltage(voltage);
    return excitation;
}

// OperatingPointExcitation Inputs::reflect_waveforms(OperatingPointExcitation excitation, double ratio) {
//     OperatingPointExcitation reflectedExcitation;
//     if (excitation.get_current()) {
//         auto reflectedCurrent = reflect_waveform(excitation.get_current().value(), ratio);
//         reflectedExcitation.set_current(reflectedCurrent);
//     }
//     if (excitation.get_voltage()) {
//         auto reflectedVoltage = reflect_waveform(excitation.get_voltage().value(), ratio);
//         reflectedExcitation.set_voltage(reflectedVoltage);
//     }
// }

SignalDescriptor Inputs::reflect_waveform(SignalDescriptor primarySignalDescriptor,
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

SignalDescriptor Inputs::reflect_waveform(SignalDescriptor signal,
                                                 double ratio,
                                                 WaveformLabel label) {
    // ABT #1670: every label is handled explicitly. An ideal transformer scales the
    // instantaneous winding voltage by 1/n and the load-referred current by n without
    // changing its shape, so every shape that is not a flyback current is reflected by
    // plain scaling. Only the flyback labels below change shape. A unipolar pulse is
    // not a flyback: a current transformer's or forward's secondary conducts DURING the
    // primary pulse, i2 = n i1 point by point (ABT #1718 review); the flyback's
    // off-time transfer has its own FLYBACK_* labels. Before,
    // any label without a case (RECTANGULAR among them) fell through a `default:`, so
    // a label nobody had thought about silently took somebody else's branch.
    switch (label) {
        case WaveformLabel::CUSTOM:
        case WaveformLabel::SINUSOIDAL:
        case WaveformLabel::TRIANGULAR:
        case WaveformLabel::TRIANGULAR_WITH_DEADTIME:
        case WaveformLabel::RECTANGULAR:
        case WaveformLabel::RECTANGULAR_DCM:
        case WaveformLabel::RECTANGULAR_WITH_DEADTIME:
        case WaveformLabel::SECONDARY_RECTANGULAR:
        case WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME:
        case WaveformLabel::BIPOLAR_RECTANGULAR:
        case WaveformLabel::BIPOLAR_TRIANGULAR:
        case WaveformLabel::UNIPOLAR_TRIANGULAR:
        case WaveformLabel::UNIPOLAR_RECTANGULAR:
            return reflect_waveform(signal, ratio);
        case WaveformLabel::FLYBACK_PRIMARY:
        case WaveformLabel::FLYBACK_SECONDARY:
            break;
        case WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME:
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "reflect_waveform: a flybackSecondaryWithDeadtime waveform has no reflection: the primary"
                " current of a discontinuous flyback depends on the on-time, which this waveform does not"
                " carry. Supply the other winding's excitation directly.");
    }

    if (!signal.get_waveform() || !signal.get_waveform()->get_time()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA,
            "reflect_waveform: reflecting a labelled waveform needs the waveform with its time axis");
    }
    // ABT #1670: the shape parameters come from the waveform being reflected, never from
    // the signal's processed block. The label is (by the bindings) guessed from the
    // waveform, and the processed block a client sends along can describe a different
    // waveform altogether: the web's default voltage carried processed {peakToPeak 100,
    // offset 0} over data -20.5/70.5, and mixing the data's label with the block's
    // numbers is how a 1:1 reflection came out as 0/-100 V.
    auto waveformToProcess = signal.get_waveform().value();
    if (is_waveform_sampled(waveformToProcess)) {
        waveformToProcess = compress_waveform(waveformToProcess);
    }
    // The label is checked on, and read off, the waveform as given, not on its compressed copy (ABT #1878).
    ProcessedWaveform processed = calculate_basic_processed_data(waveformToProcess, signal.get_waveform().value());
    Waveform newWaveform;

    double period = signal.get_waveform()->get_time()->back() - signal.get_waveform()->get_time()->front();
    double frequency = 1.0 / period;
    double peakToPeak = processed.get_peak_to_peak().value() * ratio;
    double offset = processed.get_offset() * ratio;

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
        default:
            // Unreachable for the labels dispatched above; a label added to the schema later
            // must be given its own reflection, not borrow another label's.
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "reflect_waveform: no reflection is defined for waveform label " +
                std::to_string(static_cast<int>(label)));
    }

    SignalDescriptor newSignal;
    newSignal.set_waveform(newWaveform);
    return newSignal;
}

// The excitation of the other winding of an ideal two-winding transformer: the voltage scaled
// by voltageRatio (Faraday's law, point by point, whatever the shape), the current reflected by
// its label with currentRatio, and both re-processed. Shared by calculate_reflected_secondary and
// calculate_reflected_primary so the two directions cannot drift apart.
static OperatingPointExcitation reflect_excitation(const OperatingPointExcitation& excitation,
                                                   double voltageRatio,
                                                   double currentRatio,
                                                   std::optional<std::string> name,
                                                   const std::string& functionName) {
    if (!excitation.get_voltage() || !excitation.get_voltage()->get_waveform() ||
        !excitation.get_current() || !excitation.get_current()->get_waveform()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA,
            functionName + ": the excitation needs both a voltage and a current waveform to reflect");
    }
    OperatingPointExcitation excitationOfThisWinding(excitation);
    // ABT #1670: the copy carried the other winding's name, so the secondary was presented as
    // "Primary winding excitation". It is named after its own winding when that is known,
    // and carries no name otherwise rather than the other winding's.
    excitationOfThisWinding.set_name(name);
    auto currentSignalDescriptorProcessed = Inputs::calculate_basic_processed_data(excitation.get_current().value().get_waveform().value());

    // ABT #1670: the voltage is reflected by Faraday's law alone. Every winding links the
    // same flux, so v2(t) = v1(t) / n whatever the shape; no label changes that. Routing it
    // through the label-dispatched reflection turned the web's default -20.5/70.5 V
    // rectangle into 0/-100 V, a 50 V DC level across a winding.
    auto voltageSignalDescriptor = Inputs::reflect_waveform(excitation.get_voltage().value(), voltageRatio);
    auto currentSignalDescriptor = Inputs::reflect_waveform(excitation.get_current().value(), currentRatio, currentSignalDescriptorProcessed.get_label());

    auto voltageSampledWaveform = Inputs::calculate_sampled_waveform(voltageSignalDescriptor.get_waveform().value(), excitationOfThisWinding.get_frequency());
    voltageSignalDescriptor.set_harmonics(Inputs::calculate_harmonics_data(voltageSignalDescriptor.get_waveform().value(), voltageSampledWaveform, excitationOfThisWinding.get_frequency()));
    voltageSignalDescriptor.set_processed(Inputs::calculate_processed_data(voltageSignalDescriptor, voltageSampledWaveform, true));

    auto currentSampledWaveform = Inputs::calculate_sampled_waveform(currentSignalDescriptor.get_waveform().value(), excitationOfThisWinding.get_frequency());
    currentSignalDescriptor.set_harmonics(Inputs::calculate_harmonics_data(currentSignalDescriptor.get_waveform().value(), currentSampledWaveform, excitationOfThisWinding.get_frequency()));
    currentSignalDescriptor.set_processed(Inputs::calculate_processed_data(currentSignalDescriptor, currentSampledWaveform, true));

    excitationOfThisWinding.set_voltage(voltageSignalDescriptor);
    excitationOfThisWinding.set_current(currentSignalDescriptor);

    return excitationOfThisWinding;
}

OperatingPointExcitation Inputs::calculate_reflected_secondary(OperatingPointExcitation primaryExcitation, double turnRatio, std::optional<std::string> secondaryName){
    return reflect_excitation(primaryExcitation, 1.0 / turnRatio, turnRatio, secondaryName, "calculate_reflected_secondary");
}

// ABT #1698: the web reflected the primary in its own binding, a copy of the old secondary
// logic with plain current scaling and the aliased harmonics; the reflection lives here now.
OperatingPointExcitation Inputs::calculate_reflected_primary(OperatingPointExcitation secondaryExcitation, double turnRatio, std::optional<std::string> primaryName){
    return reflect_excitation(secondaryExcitation, turnRatio, 1.0 / turnRatio, primaryName, "calculate_reflected_primary");
}

std::pair<bool, std::string> Inputs::check_integrity() {
    auto operatingPoints = get_mutable_operating_points();
    auto turnsRatios = get_design_requirements().get_turns_ratios();
    // The required magnetizing inductance is needed by ONE branch below: deriving a
    // magnetizing current for an excitation that carries a voltage but no current.
    // Resolving it eagerly here threw for every input that never reaches that branch --
    // in particular a MAS file carrying only a magnetic, whose operatingPoints array is
    // empty, so the loop body never runs. from_file computes the inductance from the
    // geometry and hands it to this constructor precisely so such a file loads, but that
    // value is only applied later in process(), leaving this line to reject the file
    // first with "DimensionWithTolerance has neither nominal, minimum nor maximum set" --
    // a complaint about a requirement nothing had asked for. Resolve where it is used, so
    // a genuinely missing inductance still throws, but only when it is genuinely needed.
    std::optional<double> resolvedMagnetizingInductance;
    auto magnetizing_inductance = [&]() {
        if (!resolvedMagnetizingInductance) {
            resolvedMagnetizingInductance =
                resolve_dimensional_values(get_design_requirements().get_magnetizing_inductance());
        }
        return resolvedMagnetizingInductance.value();
    };
    std::pair<bool, std::string> result;
    result.first = true;
    result.second = "";

    for (auto& operatingPoint : operatingPoints) {
        if (operatingPoint.get_excitations_per_winding().size() == 0) {
            result.first = false;
            throw std::invalid_argument("Missing excitation for primary");
        }
    }

    std::vector<double> turnsRatiosValues;
    for (size_t turnsRatioIndex = 0; turnsRatioIndex < turnsRatios.size(); ++turnsRatioIndex) {
        double turnsRatio = resolve_dimensional_values(turnsRatios[turnsRatioIndex]);
        turnsRatiosValues.push_back(turnsRatio);
    }

    for (size_t i = 0; i < operatingPoints.size(); ++i) {
        std::vector<OperatingPointExcitation> processedExcitationsPerWinding;
        for (auto& excitation : operatingPoints[i].get_mutable_excitations_per_winding()) {
            // Here we processed this excitation voltage
            if (excitation.get_voltage()) {
                auto voltageExcitation = excitation.get_voltage().value();
                voltageExcitation = standardize_waveform(voltageExcitation, excitation.get_frequency());
                excitation.set_voltage(voltageExcitation);
            }
            // Here we processed this excitation current
            if (excitation.get_current()) {
                auto currentExcitation = excitation.get_current().value();
                currentExcitation = standardize_waveform(currentExcitation, excitation.get_frequency());
                excitation.set_current(currentExcitation);
            }
            else {
                // ABT #825: with no current the current is derived from the voltage — which was
                // dereferenced unconditionally, twice. An excitation carrying neither (or a voltage
                // with no waveform) therefore died on an empty optional, and "bad optional access"
                // named neither the operating point, the winding, nor what was actually needed.
                if (!excitation.get_voltage() || !excitation.get_voltage()->get_waveform()) {
                    std::string operatingPointLabel = "operating point " + std::to_string(i) +
                        (operatingPoints[i].get_name() ? " ('" + operatingPoints[i].get_name().value() + "')" : "");
                    throw std::invalid_argument(
                        operatingPointLabel + ", winding " + std::to_string(processedExcitationsPerWinding.size()) +
                        ": the excitation has no current, so one has to be derived from the voltage — and it has no"
                        " voltage waveform either. Give the winding a current, or a voltage to derive it from.");
                }
                auto voltageWaveform = excitation.get_voltage()->get_waveform().value();
                auto sampledWaveform = calculate_sampled_waveform(voltageWaveform, excitation.get_frequency());
                // This winding has no current, so no current gives its magnetizing DC: Faraday alone.
                excitation.set_current(integrate_magnetizing_current(excitation, sampledWaveform, magnetizing_inductance(), true, 0.0));
            }
            processedExcitationsPerWinding.push_back(excitation);
        }
        operatingPoints[i].set_excitations_per_winding(processedExcitationsPerWinding);
    }

    for (size_t i = 0; i < operatingPoints.size(); ++i) {

        // ABT #825: the secondary is synthesised by reflecting the primary, and the primary's
        // voltage and current were both dereferenced unconditionally. An operating point that
        // gives a current and no voltage — a perfectly ordinary way to describe an inductor's
        // winding — therefore died on an empty optional, with "bad optional access" naming
        // neither the operating point, nor the winding, nor the field. Reflect whatever the
        // primary actually has, and when it has nothing to reflect, say so.
        if (turnsRatios.size() > operatingPoints[i].get_excitations_per_winding().size() - 1) {
            std::string operatingPointLabel = "operating point " + std::to_string(i) +
                (operatingPoints[i].get_name() ? " ('" + operatingPoints[i].get_name().value() + "')" : "");

            if (turnsRatios.size() == 1 && operatingPoints[i].get_excitations_per_winding().size() == 1) {
                // We are missing excitation only for secondary
                for (size_t turnsRatioIndex = 0; turnsRatioIndex < turnsRatios.size(); ++turnsRatioIndex) {
                    if (turnsRatioIndex >= operatingPoints[i].get_excitations_per_winding().size() - 1) {
                        double turnsRatio = resolve_dimensional_values(turnsRatios[turnsRatioIndex]);
                        auto excitationOfPrimaryWinding = operatingPoints[i].get_excitations_per_winding()[0];

                        if (!excitationOfPrimaryWinding.get_voltage() && !excitationOfPrimaryWinding.get_current()) {
                            throw std::invalid_argument(
                                operatingPointLabel + ": the design requirements declare " +
                                std::to_string(turnsRatios.size() + 1) + " windings but only one excitation was"
                                " given, so the secondary has to be reflected from the primary — and winding 0"
                                " has neither a voltage nor a current to reflect. Give winding 0 a voltage or a"
                                " current, or supply the secondary excitation directly.");
                        }

                        OperatingPointExcitation excitationOfThisWinding(excitationOfPrimaryWinding);
                        // ABT #1670: the copy is not the primary; it must not carry its name.
                        excitationOfThisWinding.set_name(std::nullopt);

                        // Reflect only what the primary actually carries. Inventing the other half
                        // would be making up an operating point nobody asked for.
                        if (excitationOfPrimaryWinding.get_voltage()) {
                            excitationOfThisWinding.set_voltage(
                                reflect_waveform(excitationOfPrimaryWinding.get_voltage().value(), 1 / turnsRatio));
                        }
                        if (excitationOfPrimaryWinding.get_current()) {
                            excitationOfThisWinding.set_current(
                                reflect_waveform(excitationOfPrimaryWinding.get_current().value(), turnsRatio));
                        }
                        operatingPoints[i].get_mutable_excitations_per_winding().push_back(excitationOfThisWinding);
                    }
                }
                result.second = "Had to create the excitations of some windings based on primary";
            }
            else {
                throw std::invalid_argument(
                    operatingPointLabel + ": the design requirements declare " +
                    std::to_string(turnsRatios.size() + 1) + " windings but only " +
                    std::to_string(operatingPoints[i].get_excitations_per_winding().size()) +
                    " excitation(s) were given. Only a single missing secondary can be reflected from the"
                    " primary; supply an excitation per winding.");
            }
        }
    }

    set_operating_points(operatingPoints);

    return result;
}

ProcessedWaveform Inputs::calculate_processed_data(SignalDescriptor excitation,
                                            Waveform sampledWaveform,
                                            bool includeAdvancedData,
                                            std::optional<ProcessedWaveform> processed) {
    return WaveformProcessor::calculate_processed_data(excitation, sampledWaveform, includeAdvancedData, processed,
                                                       settings.get_inputs_number_points_sampled_waveforms());
}

ProcessedWaveform Inputs::calculate_basic_processed_data(Waveform waveform, std::optional<Waveform> givenWaveform) {
    // Delegates to WaveformProcessor, which additionally fixed part of ABT
    // #602 this copy never had: the duty-cycle measurement for a CUSTOM label
    // (i.e. real imported/simulated data) needs the UNCOMPRESSED waveform —
    // this copy always passed the compressed one, so real SPICE data (which
    // compress_waveform collapses to 50-60 points instead of an analytical
    // 3/4/5-point shape) had its duty measured against a thinned-out signal.
    return WaveformProcessor::calculate_basic_processed_data(waveform, settings.get_inputs_number_points_sampled_waveforms(), givenWaveform);
}

ProcessedWaveform Inputs::calculate_processed_data(Waveform waveform,
                                            double frequency,
                                            bool includeAdvancedData,
                                            std::optional<ProcessedWaveform> processed) {
    return WaveformProcessor::calculate_processed_data(waveform, frequency, includeAdvancedData, processed,
                                                       settings.get_inputs_trim_harmonics(),
                                                       settings.get_harmonic_amplitude_threshold(),
                                                       settings.get_inputs_number_points_sampled_waveforms());
}

ProcessedWaveform Inputs::calculate_processed_data(Harmonics harmonics,
                                            Waveform waveform,
                                            bool includeAdvancedData,
                                            std::optional<ProcessedWaveform> processed) {
    return WaveformProcessor::calculate_processed_data(harmonics, waveform, includeAdvancedData, processed,
                                                       settings.get_inputs_number_points_sampled_waveforms());
}

Harmonics Inputs::calculate_harmonics_data(Waveform waveform, double frequency) {
    return WaveformProcessor::calculate_harmonics_data(waveform, frequency,
                                                       settings.get_inputs_trim_harmonics(),
                                                       settings.get_harmonic_amplitude_threshold(),
                                                       settings.get_inputs_number_points_sampled_waveforms());
}

Harmonics Inputs::calculate_harmonics_data(const Waveform& waveform, const Waveform& sampledWaveform, double frequency) {
    return calculate_harmonics_data(waveform.get_time() ? waveform : sampledWaveform, frequency);
}

OperatingPointExcitation Inputs::prune_harmonics(OperatingPointExcitation excitation, double windingLossesHarmonicAmplitudeThreshold, std::optional<size_t> mainHarmonicIndex) {
    if (excitation.get_current()) {
        excitation.set_current(prune_harmonics(excitation.get_current().value(), windingLossesHarmonicAmplitudeThreshold, mainHarmonicIndex));
    }
    if (excitation.get_voltage()) {
        excitation.set_voltage(prune_harmonics(excitation.get_voltage().value(), windingLossesHarmonicAmplitudeThreshold, mainHarmonicIndex));
    }
    if (excitation.get_magnetizing_current()) {
        excitation.set_magnetizing_current(prune_harmonics(excitation.get_magnetizing_current().value(), windingLossesHarmonicAmplitudeThreshold, mainHarmonicIndex));
    }

    return excitation;
}

SignalDescriptor Inputs::prune_harmonics(SignalDescriptor signalDescriptor, double windingLossesHarmonicAmplitudeThreshold, std::optional<size_t> mainHarmonicIndex) {
    if (!signalDescriptor.get_harmonics()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Signal has no harmonics to prune");
    }
    auto unprunedHarmonics = signalDescriptor.get_harmonics().value();
    Harmonics prunedHarmonics;
    auto harmonicsIndexesToMaintain = get_main_harmonic_indexes(unprunedHarmonics, windingLossesHarmonicAmplitudeThreshold, mainHarmonicIndex);
    prunedHarmonics.get_mutable_amplitudes().push_back(unprunedHarmonics.get_amplitudes()[0]);
    prunedHarmonics.get_mutable_frequencies().push_back(unprunedHarmonics.get_frequencies()[0]);

    for (auto harmonicIndex : harmonicsIndexesToMaintain) {
        prunedHarmonics.get_mutable_amplitudes().push_back(unprunedHarmonics.get_amplitudes()[harmonicIndex]);
        prunedHarmonics.get_mutable_frequencies().push_back(unprunedHarmonics.get_frequencies()[harmonicIndex]);
    }

    signalDescriptor.set_harmonics(prunedHarmonics);

    return signalDescriptor;    
}

OperatingPoint Inputs::prune_harmonics(OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold, std::optional<size_t> mainHarmonicIndex) {
    bool allExcitationsHaveCurrent = true;
    bool allExcitationsHaveMagnetizingCurrent = true;
    bool allExcitationsHaveVoltage = true;

    for (auto excitation : operatingPoint.get_excitations_per_winding()) {
        // A signal only takes part in pruning if it actually carries harmonics;
        // the branches below dereference get_harmonics().value() (ABT #599: a
        // SignalDescriptor present without harmonics threw bad_optional_access).
        if (!excitation.get_current() || !excitation.get_current()->get_harmonics()) {
            allExcitationsHaveCurrent = false;
        }
        if (!excitation.get_magnetizing_current() || !excitation.get_magnetizing_current()->get_harmonics()) {
            allExcitationsHaveMagnetizingCurrent = false;
        }
        if (!excitation.get_voltage() || !excitation.get_voltage()->get_harmonics()) {
            allExcitationsHaveVoltage = false;
        }
    }
    std::vector<size_t> currentHarmonicIndexesToMaintain;
    std::vector<double> currentHarmonicFrequenciesToMaintain;
    std::vector<size_t> magnetizingCurrentHarmonicIndexesToMaintain;
    std::vector<double> magnetizingCurrentHarmonicFrequenciesToMaintain;
    std::vector<size_t> voltageHarmonicIndexesToMaintain;
    std::vector<double> voltageHarmonicFrequenciesToMaintain;

    for (auto excitation : operatingPoint.get_excitations_per_winding()) {
        // Current
        if (allExcitationsHaveCurrent) {
            auto unprunedHarmonics = excitation.get_current()->get_harmonics().value();
            auto harmonicsIndexesToMaintain = get_main_harmonic_indexes(unprunedHarmonics, windingLossesHarmonicAmplitudeThreshold, mainHarmonicIndex);
            for (auto harmonicIndex : harmonicsIndexesToMaintain) {
                if (std::find(currentHarmonicIndexesToMaintain.begin(), currentHarmonicIndexesToMaintain.end(), harmonicIndex) == currentHarmonicIndexesToMaintain.end()) {
                    currentHarmonicIndexesToMaintain.push_back(harmonicIndex);
                    currentHarmonicFrequenciesToMaintain.push_back(unprunedHarmonics.get_frequencies()[harmonicIndex]);
                }
            }
        }

        // Magnetizing current
        if (allExcitationsHaveMagnetizingCurrent) {
            auto unprunedHarmonics = excitation.get_magnetizing_current()->get_harmonics().value();
            auto harmonicsIndexesToMaintain = get_main_harmonic_indexes(unprunedHarmonics, windingLossesHarmonicAmplitudeThreshold, mainHarmonicIndex);
            for (auto harmonicIndex : harmonicsIndexesToMaintain) {
                if (std::find(magnetizingCurrentHarmonicIndexesToMaintain.begin(), magnetizingCurrentHarmonicIndexesToMaintain.end(), harmonicIndex) == magnetizingCurrentHarmonicIndexesToMaintain.end()) {
                    magnetizingCurrentHarmonicIndexesToMaintain.push_back(harmonicIndex);
                    magnetizingCurrentHarmonicFrequenciesToMaintain.push_back(unprunedHarmonics.get_frequencies()[harmonicIndex]);
                }
            }
        }

        // Voltage
        if (allExcitationsHaveVoltage) {
            auto unprunedHarmonics = excitation.get_voltage()->get_harmonics().value();
            auto harmonicsIndexesToMaintain = get_main_harmonic_indexes(unprunedHarmonics, windingLossesHarmonicAmplitudeThreshold, mainHarmonicIndex);
            for (auto harmonicIndex : harmonicsIndexesToMaintain) {
                if (std::find(voltageHarmonicIndexesToMaintain.begin(), voltageHarmonicIndexesToMaintain.end(), harmonicIndex) == voltageHarmonicIndexesToMaintain.end()) {
                    voltageHarmonicIndexesToMaintain.push_back(harmonicIndex);
                    voltageHarmonicFrequenciesToMaintain.push_back(unprunedHarmonics.get_frequencies()[harmonicIndex]);
                }
            }
        }
    }

    std::vector<OperatingPointExcitation> newExcitationsPerWinding;
    for (auto excitation : operatingPoint.get_excitations_per_winding()) {
        // Start from a full copy of the original excitation so every field that
        // pruning does not touch (frequency, magneticFluxDensity,
        // magneticFieldStrength, ...) survives. Rebuilding field-by-field
        // silently dropped `frequency` (left at its default of 0), which zeroed
        // the frequency of every analytical operating point produced by
        // Topology::process() and broke all downstream frequency-derived math.
        OperatingPointExcitation newExcitation = excitation;
        // Current
        if (allExcitationsHaveCurrent) {
            auto current = excitation.get_current().value();
            auto unprunedHarmonics = current.get_harmonics().value();
            Harmonics prunedHarmonics;
            prunedHarmonics.get_mutable_amplitudes().push_back(unprunedHarmonics.get_amplitudes()[0]);
            prunedHarmonics.get_mutable_frequencies().push_back(unprunedHarmonics.get_frequencies()[0]);
            // Sort indexes and frequencies together to maintain alignment
            {
                std::vector<std::pair<size_t, double>> indexFreqPairs;
                for (size_t i = 0; i < currentHarmonicIndexesToMaintain.size(); ++i) {
                    indexFreqPairs.push_back({currentHarmonicIndexesToMaintain[i], currentHarmonicFrequenciesToMaintain[i]});
                }
                std::sort(indexFreqPairs.begin(), indexFreqPairs.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
                currentHarmonicIndexesToMaintain.clear();
                currentHarmonicFrequenciesToMaintain.clear();
                for (const auto& pair : indexFreqPairs) {
                    currentHarmonicIndexesToMaintain.push_back(pair.first);
                    currentHarmonicFrequenciesToMaintain.push_back(pair.second);
                }
            }
            for (size_t i = 0; i < currentHarmonicIndexesToMaintain.size(); ++i) {
                auto harmonicIndex = currentHarmonicIndexesToMaintain[i];
                if (harmonicIndex < unprunedHarmonics.get_amplitudes().size()) {
                    prunedHarmonics.get_mutable_amplitudes().push_back(unprunedHarmonics.get_amplitudes()[harmonicIndex]);
                    prunedHarmonics.get_mutable_frequencies().push_back(unprunedHarmonics.get_frequencies()[harmonicIndex]);
                }
                else {
                    prunedHarmonics.get_mutable_amplitudes().push_back(0);
                    prunedHarmonics.get_mutable_frequencies().push_back(currentHarmonicFrequenciesToMaintain[i]);
                }
            }
            current.set_harmonics(prunedHarmonics);
            newExcitation.set_current(current);
        }
        else {
            newExcitation.set_current(excitation.get_current());
        }

        // Magnetizing current
        if (allExcitationsHaveMagnetizingCurrent) {
            auto magnetizingCurrent = excitation.get_magnetizing_current().value();
            auto unprunedHarmonics = magnetizingCurrent.get_harmonics().value();
            Harmonics prunedHarmonics;
            prunedHarmonics.get_mutable_amplitudes().push_back(unprunedHarmonics.get_amplitudes()[0]);
            prunedHarmonics.get_mutable_frequencies().push_back(unprunedHarmonics.get_frequencies()[0]);
            // Sort indexes and frequencies together to maintain alignment
            {
                std::vector<std::pair<size_t, double>> indexFreqPairs;
                for (size_t i = 0; i < magnetizingCurrentHarmonicIndexesToMaintain.size(); ++i) {
                    indexFreqPairs.push_back({magnetizingCurrentHarmonicIndexesToMaintain[i], magnetizingCurrentHarmonicFrequenciesToMaintain[i]});
                }
                std::sort(indexFreqPairs.begin(), indexFreqPairs.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
                magnetizingCurrentHarmonicIndexesToMaintain.clear();
                magnetizingCurrentHarmonicFrequenciesToMaintain.clear();
                for (const auto& pair : indexFreqPairs) {
                    magnetizingCurrentHarmonicIndexesToMaintain.push_back(pair.first);
                    magnetizingCurrentHarmonicFrequenciesToMaintain.push_back(pair.second);
                }
            }
            for (size_t i = 0; i < magnetizingCurrentHarmonicIndexesToMaintain.size(); ++i) {
                auto harmonicIndex = magnetizingCurrentHarmonicIndexesToMaintain[i];
                if (harmonicIndex < unprunedHarmonics.get_amplitudes().size()) {
                    prunedHarmonics.get_mutable_amplitudes().push_back(unprunedHarmonics.get_amplitudes()[harmonicIndex]);
                    prunedHarmonics.get_mutable_frequencies().push_back(unprunedHarmonics.get_frequencies()[harmonicIndex]);
                }
                else {
                    prunedHarmonics.get_mutable_amplitudes().push_back(0);
                    prunedHarmonics.get_mutable_frequencies().push_back(magnetizingCurrentHarmonicFrequenciesToMaintain[i]);
                }
            }
            magnetizingCurrent.set_harmonics(prunedHarmonics);
            newExcitation.set_magnetizing_current(magnetizingCurrent);
        }
        else {
            newExcitation.set_magnetizing_current(excitation.get_magnetizing_current());
        }


        // Voltage
        if (allExcitationsHaveVoltage) {
            auto voltage = excitation.get_voltage().value();
            auto unprunedHarmonics = voltage.get_harmonics().value();
            Harmonics prunedHarmonics;
            prunedHarmonics.get_mutable_amplitudes().push_back(unprunedHarmonics.get_amplitudes()[0]);
            prunedHarmonics.get_mutable_frequencies().push_back(unprunedHarmonics.get_frequencies()[0]);
            // Sort indexes and frequencies together to maintain alignment
            {
                std::vector<std::pair<size_t, double>> indexFreqPairs;
                for (size_t i = 0; i < voltageHarmonicIndexesToMaintain.size(); ++i) {
                    indexFreqPairs.push_back({voltageHarmonicIndexesToMaintain[i], voltageHarmonicFrequenciesToMaintain[i]});
                }
                std::sort(indexFreqPairs.begin(), indexFreqPairs.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
                voltageHarmonicIndexesToMaintain.clear();
                voltageHarmonicFrequenciesToMaintain.clear();
                for (const auto& pair : indexFreqPairs) {
                    voltageHarmonicIndexesToMaintain.push_back(pair.first);
                    voltageHarmonicFrequenciesToMaintain.push_back(pair.second);
                }
            }
            for (size_t i = 0; i < voltageHarmonicIndexesToMaintain.size(); ++i) {
                auto harmonicIndex = voltageHarmonicIndexesToMaintain[i];
                if (harmonicIndex < unprunedHarmonics.get_amplitudes().size()) {
                    prunedHarmonics.get_mutable_amplitudes().push_back(unprunedHarmonics.get_amplitudes()[harmonicIndex]);
                    prunedHarmonics.get_mutable_frequencies().push_back(unprunedHarmonics.get_frequencies()[harmonicIndex]);
                }
                else {
                    prunedHarmonics.get_mutable_amplitudes().push_back(0);
                    prunedHarmonics.get_mutable_frequencies().push_back(voltageHarmonicFrequenciesToMaintain[i]);
                }
            }
            voltage.set_harmonics(prunedHarmonics);
            newExcitation.set_voltage(voltage);
        }
        else {
            newExcitation.set_voltage(excitation.get_voltage());
        }

        // FIXED: BUG-03 — Push pruned newExcitation, not original excitation
        newExcitationsPerWinding.push_back(newExcitation);
    }
    operatingPoint.set_excitations_per_winding(newExcitationsPerWinding);

    return operatingPoint;    
}


Waveform Inputs::compress_waveform(const Waveform& waveform, std::optional<double> frequency) {
    return WaveformProcessor::compress_waveform(waveform, frequency);
}

Inputs::VoltSecondBalance Inputs::calculate_volt_second_balance(const Waveform& voltageSampledWaveform, const std::string& context) {
    // ABT #1585: the magnetizing current is the integral of the voltage over the inductance. Only
    // the integration constant is removed (calculate_integral_waveform), so a voltage whose
    // volt-seconds do not balance over the period ramps the "current" by its imbalance instead of
    // describing the flux swing: Kirchhoff's line-cycle PFC inductor voltages gave 1635 A against
    // 12.3 A of real winding current, and every core then failed the energy-stored gate.
    if (!voltageSampledWaveform.get_time() || voltageSampledWaveform.get_data().size() < 2) {
        throw InvalidInputException(ErrorCode::MISSING_DATA,
            context + ": the voltage needs a sampled waveform with a time axis to check its volt-second balance");
    }
    VoltSecondBalance balance;
    balance.meanVoltage = calculate_waveform_average(voltageSampledWaveform);
    const auto time = voltageSampledWaveform.get_time().value();
    balance.period = time.back() - time.front();
    balance.netVoltSeconds = balance.meanVoltage * balance.period;
    auto acFluxLinkage = calculate_integral_waveform(sum_waveform(voltageSampledWaveform, -balance.meanVoltage), true);
    const auto& acData = acFluxLinkage.get_data();
    balance.acSwing = *std::max_element(acData.begin(), acData.end()) - *std::min_element(acData.begin(), acData.end());
    if (balance.acSwing > 0) {
        balance.imbalance = std::fabs(balance.netVoltSeconds) / balance.acSwing;
    }
    else {
        // A constant voltage: any non-zero value is pure imbalance; zero is no voltage at all.
        balance.imbalance = balance.meanVoltage != 0 ? std::numeric_limits<double>::infinity() : 0;
    }
    return balance;
}

std::optional<std::string> Inputs::describe_volt_second_imbalance(const Waveform& voltageSampledWaveform, const std::string& context) {
    auto balance = calculate_volt_second_balance(voltageSampledWaveform, context);
    if (balance.acSwing <= 0) {
        if (balance.meanVoltage != 0) {
            return context + ": the voltage is a constant " + format_quantity(balance.meanVoltage, "V") +
                ". A constant voltage has no periodic steady state, so its integral is not a magnetizing current";
        }
        return std::nullopt;
    }
    if (balance.imbalance > maximumVoltSecondImbalance) {
        return context + ": the voltage is not volt-second balanced: its mean of " + format_quantity(balance.meanVoltage, "V") +
            " drifts the flux linkage by " + format_quantity(balance.netVoltSeconds, "V*s") + " every " +
            format_quantity(balance.period, "s") + " period, " + format_significant(balance.imbalance, 2) +
            " times its AC swing of " + format_quantity(balance.acSwing, "V*s") + " (the limit is " +
            format_significant(maximumVoltSecondImbalance, 2) + "). A magnetizing current integrated from it would "
            "describe the drift: give a balanced steady-state voltage, or the magnetizing current itself";
    }
    return std::nullopt;
}

void Inputs::check_volt_second_balance(const Waveform& voltageSampledWaveform, const std::string& context) {
    auto imbalance = describe_volt_second_imbalance(voltageSampledWaveform, context);
    if (imbalance) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, imbalance.value());
    }
}

namespace {

// The left and right limits of a piecewise-linear waveform at `time`, inside its time axis: a vertical
// edge (repeated time) has a value on each side, any other instant one value.
std::pair<double, double> waveform_limits_at(const std::vector<double>& time, const std::vector<double>& data, double instant, double tolerance) {
    auto first = std::lower_bound(time.begin(), time.end(), instant - tolerance);
    auto last = std::upper_bound(time.begin(), time.end(), instant + tolerance);
    if (first != last) {
        return {data[first - time.begin()], data[(last - time.begin()) - 1]};
    }
    if (first == time.begin() || first == time.end()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "remove_resistive_drop: instant " + format_quantity(instant, "s") +
                                    " is outside the waveform's time axis");
    }
    size_t after = first - time.begin();
    size_t before = after - 1;
    double fraction = (instant - time[before]) / (time[after] - time[before]);
    double value = data[before] + fraction * (data[after] - data[before]);
    return {value, value};
}

// The same limits for a periodic waveform at any instant: `instant` is mapped into the waveform's own
// period [t0, t0 + period], and the stretch between its last knot and t0 + period (an axis that stops
// one sample short of the period) runs linearly back to its first value. At the period boundary the
// left limit is the end of the period and the right limit its start.
std::pair<double, double> periodic_waveform_limits_at(const std::vector<double>& time, const std::vector<double>& data, double instant,
                                                      double period, double tolerance) {
    auto limits_at_shift = [&](double shift) -> std::pair<double, double> {
        double mapped = time.front() + shift;
        if (mapped <= time.back() + tolerance) {
            return waveform_limits_at(time, data, std::min(mapped, time.back()), tolerance);
        }
        double fraction = (mapped - time.back()) / (time.front() + period - time.back());
        double value = data.back() + fraction * (data.front() - data.back());
        return {value, value};
    };
    double shift = std::fmod(instant - time.front(), period);
    if (shift < 0) {
        shift += period;
    }
    if (period - shift <= tolerance) {
        shift = 0;
    }
    if (shift <= tolerance) {
        auto start = limits_at_shift(0);
        auto end = limits_at_shift(period);
        return {end.first, start.second};
    }
    return limits_at_shift(shift);
}

} // namespace

OperatingPointExcitation Inputs::remove_resistive_drop(const OperatingPointExcitation& excitation, double dcResistance) {
    std::string name = "Excitation '" + excitation.get_name().value_or("unnamed") + "'";
    if (!std::isfinite(dcResistance) || dcResistance < 0) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, name + ": a winding's DC resistance must be a finite value >= 0, not " +
                                    format_quantity(dcResistance, "Ohm"));
    }
    if (!excitation.get_voltage()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, name + ": no voltage to remove the resistive drop from");
    }
    if (!excitation.get_current()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, name + ": no current, so its resistive drop is not known");
    }
    double frequency = excitation.get_frequency();
    if (!(frequency > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, name + ": the frequency must be positive to read its current periodically");
    }
    double period = 1.0 / frequency;
    auto voltage = standardize_waveform(excitation.get_voltage().value(), frequency).get_waveform().value();
    auto current = standardize_waveform(excitation.get_current().value(), frequency).get_waveform().value();
    const auto voltageTime = voltage.get_time().value();
    const auto& voltageData = voltage.get_data();
    const auto currentTime = current.get_time().value();
    const auto& currentData = current.get_data();
    if (voltageTime.size() < 2 || currentTime.size() < 2 || voltageTime.size() != voltageData.size() || currentTime.size() != currentData.size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, name + ": the voltage and the current need at least two points with one time each");
    }
    double tolerance = period * 1e-9;
    if (currentTime.back() - currentTime.front() > period + tolerance * 1e3) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, name + ": the current spans " + format_quantity(currentTime.back() - currentTime.front(), "s") +
                                    ", more than one period (" + format_quantity(period, "s") + "), so it cannot be read periodically");
    }

    // The voltage must be given over one period, on a closed axis [t0, t0 + T] or on an open sampled one
    // [t0, t0 + T - dt] (WaveformProcessor::check_time_axis_spans_one_period, which throws otherwise). Both
    // are read periodically, an open axis running linearly from its last sample back to its first one at
    // t0 + T, as the sampler reads it. The corrected voltage keeps the voltage's own axis when it can: a
    // sampled voltage whose current has no knot and no edge off its samples (both sampled on one grid)
    // stays the same sampled waveform, so it is sampled to the same number of points as the voltage it
    // replaces. Otherwise it is given on the closed axis [t0, t0 + T], with the knots of both waveforms.
    const bool openVoltageAxis = WaveformProcessor::check_time_axis_spans_one_period(voltageTime, frequency);
    const double periodStart = voltageTime.front();
    const double periodEnd = periodStart + period;
    std::vector<double> instants;
    for (double instant : voltageTime) {
        // A closed axis ends at t0 + T within the period check's tolerance: that end is periodEnd.
        if (instant < periodEnd - tolerance) {
            instants.push_back(instant);
        }
    }
    const size_t numberVoltageKnots = instants.size();
    for (double instant : currentTime) {
        double shift = std::fmod(instant - periodStart, period);
        if (shift < 0) {
            shift += period;
        }
        if (period - shift > tolerance) {
            instants.push_back(periodStart + shift);
        }
    }
    std::sort(instants.begin(), instants.end());
    std::vector<double> uniqueInstants;
    for (double instant : instants) {
        if (uniqueInstants.empty() || instant - uniqueInstants.back() > tolerance) {
            uniqueInstants.push_back(instant);
        }
    }
    // Instants within tolerance of t0 + T are that period end, which the closed axis adds below.
    if (periodEnd - uniqueInstants.back() <= tolerance) {
        uniqueInstants.pop_back();
    }

    std::vector<std::pair<double, double>> limits;
    double voltagePeak = 0;
    for (double value : voltageData) {
        voltagePeak = std::max(voltagePeak, std::fabs(value));
    }
    bool singleValued = true;
    for (double instant : uniqueInstants) {
        auto [voltageLeft, voltageRight] = periodic_waveform_limits_at(voltageTime, voltageData, instant, period, tolerance);
        auto [currentLeft, currentRight] = periodic_waveform_limits_at(currentTime, currentData, instant, period, tolerance);
        double left = voltageLeft - dcResistance * currentLeft;
        double right = voltageRight - dcResistance * currentRight;
        limits.emplace_back(left, right);
        if (std::fabs(right - left) > 1e-12 * std::max({voltagePeak, std::fabs(left), std::fabs(right)})) {
            singleValued = false;
        }
    }

    std::vector<double> time;
    std::vector<double> data;
    if (openVoltageAxis && uniqueInstants.size() == numberVoltageKnots && singleValued) {
        // The voltage's own sampled axis, one value per sample.
        for (size_t index = 0; index < uniqueInstants.size(); ++index) {
            time.push_back(voltageTime[index]);
            data.push_back(limits[index].second);
        }
    }
    else {
        // The closed axis [t0, t0 + T]: the first knot starts the period, t0 + T ends it with the left
        // limit at the period boundary, and every edge inside keeps its value on both sides.
        for (size_t index = 0; index < uniqueInstants.size(); ++index) {
            auto [left, right] = limits[index];
            if (index > 0) {
                time.push_back(uniqueInstants[index]);
                data.push_back(left);
            }
            if (index == 0 || right != left) {
                time.push_back(uniqueInstants[index]);
                data.push_back(right);
            }
        }
        time.push_back(periodEnd);
        data.push_back(limits[0].first);
    }

    Waveform magnetizingVoltageWaveform;
    magnetizingVoltageWaveform.set_time(time);
    magnetizingVoltageWaveform.set_data(data);
    SignalDescriptor magnetizingVoltage;
    magnetizingVoltage.set_waveform(magnetizingVoltageWaveform);
    OperatingPointExcitation result = excitation;
    result.set_voltage(magnetizingVoltage);
    return result;
}

std::optional<std::vector<double>> Inputs::get_winding_dc_resistances(const Coil& coil, double temperature) {
    if (!coil.get_turns_description()) {
        return std::nullopt;
    }
    return WindingOhmicLosses::calculate_dc_resistance_per_winding(coil, temperature);
}

double Inputs::calculate_steady_state_mean_voltage(const OperatingPointExcitation& excitation,
                                                   const Waveform& voltageSampledWaveform,
                                                   const std::string& context) {
    // The sampled voltage must have a time axis (the same refusal as the balance check), but an
    // unbalanced voltage is not refused here: its mean is removed like any other (ABT #1855).
    calculate_volt_second_balance(voltageSampledWaveform, context);
    // The mean is taken from the voltage as given, not from its resampled copy: resampling a clean
    // rectangle moves its edges by a fraction of a sample and invents a mean of ~1 % of the swing,
    // which is not in the data and must not bend the integral.
    const Waveform& givenWaveform = (excitation.get_voltage() && excitation.get_voltage()->get_waveform() &&
                                     excitation.get_voltage()->get_waveform()->get_time() &&
                                     excitation.get_voltage()->get_waveform()->get_data().size() >= 2)
                                        ? excitation.get_voltage()->get_waveform().value()
                                        : voltageSampledWaveform;
    double meanVoltage = calculate_waveform_average(givenWaveform, excitation.get_frequency());
    if (meanVoltage == 0) {
        return 0;
    }
    const auto givenTime = givenWaveform.get_time().value();
    double period = givenTime.back() - givenTime.front();
    double meanSampledVoltage = calculate_waveform_average(voltageSampledWaveform, excitation.get_frequency());
    auto acFluxLinkage = calculate_integral_waveform(sum_waveform(voltageSampledWaveform, -meanSampledVoltage), true);
    const auto& acData = acFluxLinkage.get_data();
    double acSwing = *std::max_element(acData.begin(), acData.end()) - *std::min_element(acData.begin(), acData.end());
    double imbalance = acSwing > 0 ? std::fabs(meanVoltage * period) / acSwing : 0;
    if (imbalance > voltSecondImbalanceWarningThreshold) {
        Logger::getInstance().warning(context + ": the voltage has a mean of " + format_quantity(meanVoltage, "V") +
            ", whose net volt-seconds are " + format_significant(imbalance * 100, 2) + " % of its AC swing. MKF removed it "
            "before integrating the magnetizing current, since a winding in periodic steady state has none (it comes "
            "from a probe offset, a resistive drop or a capture that is not a whole number of periods)", "Inputs");
    }
    return meanVoltage;
}

void Inputs::check_adviser_inputs_volt_second_balance(const Inputs& inputs,
                                                      std::optional<std::vector<IsolationSide>> coilIsolationSides,
                                                      std::optional<std::string> primaryWindingName,
                                                      std::optional<CoreShapeFamily> coreShapeFamily,
                                                      std::optional<CoreType> coreType) {
    // Open-core families are priced by their own models (calculate_flux_density_for_family_model),
    // which never integrate the voltage. DRUM_RING and MOLDED take those models only for some
    // grades, so the whole family is left unchecked rather than guessing which.
    if (coreShapeFamily) {
        switch (coreShapeFamily.value()) {
            case CoreShapeFamily::DRUM:
            case CoreShapeFamily::DRUM_SEMISHIELDED:
            case CoreShapeFamily::DRUM_RING:
            case CoreShapeFamily::ROD:
            case CoreShapeFamily::MOLDED:
                return;
            default:
                break;
        }
    }

    const auto& operatingPoints = inputs.get_operating_points();
    for (size_t operatingPointIndex = 0; operatingPointIndex < operatingPoints.size(); ++operatingPointIndex) {
        const auto& operatingPoint = operatingPoints[operatingPointIndex];
        const auto& excitations = operatingPoint.get_excitations_per_winding();
        if (excitations.empty()) {
            continue;
        }
        // The loss path integrates only the primary (winding 0) excitation.
        auto excitation = get_primary_excitation(operatingPoint);
        if (!excitation.get_voltage()) {
            continue;
        }
        bool hasCurrent = excitation.get_current().has_value();

        // Mirrors the routing in MagnetizingInductance::calculate_inductance_and_magnetic_flux_density.
        // One winding in the circuit with a current: that current is the magnetizing current (the
        // loss path also requires the coil to have one winding; with the coil unknown, or with it
        // known, this check stays on the side that does not integrate).
        if (excitations.size() == 1 && hasCurrent) {
            continue;
        }
        if (excitations.size() >= 2) {
            // With no primary current, the routing below throws its own error before integrating.
            if (!hasCurrent) {
                continue;
            }
            // Common-mode choke on a toroid: its magnetizing current is the winding average. With
            // the core type unknown, every CMC-shaped operating point is skipped.
            bool coreCanBeToroidal = !coreType || coreType.value() == CoreType::TOROIDAL;
            if (coreCanBeToroidal && can_be_common_mode_choke(operatingPoint)) {
                continue;
            }
            // Multiport inductor: the coil's isolation sides decide, unless the coil has fewer
            // windings than the circuit. With no coil yet, either reading skips the check.
            std::optional<std::vector<IsolationSide>> effectiveIsolationSides = std::nullopt;
            if (coilIsolationSides && coilIsolationSides->size() >= excitations.size()) {
                effectiveIsolationSides = coilIsolationSides;
            }
            // A flyback-labelled multiport inductor has its primary voltage integrated (ABT #1467);
            // one with every winding on one side takes the ampere-turn sum and integrates nothing.
            if (!has_flyback_label(excitation)) {
                if (is_multiport_inductor(operatingPoint, effectiveIsolationSides)) {
                    continue;
                }
                if (!coilIsolationSides) {
                    auto requiredIsolationSides = inputs.get_design_requirements().get_isolation_sides();
                    if (requiredIsolationSides && is_multiport_inductor(operatingPoint, requiredIsolationSides.value())) {
                        continue;
                    }
                }
            }
        }
        // Every current label has its voltage integrated (ABT #1467): a FLYBACK_PRIMARY /
        // FLYBACK_SECONDARY / UNIPOLAR_TRIANGULAR current is no longer rebuilt from its own knots.

        // The same preparation the loss path gives the voltage before integrating it, so the
        // numbers in the message are the ones it would have reported.
        double frequency = excitation.get_frequency();
        OperatingPoint preparedOperatingPoint = operatingPoint;
        if (excitation.get_current()) {
            excitation.set_current(standardize_waveform(excitation.get_current().value(), frequency));
        }
        excitation.set_voltage(standardize_waveform(excitation.get_voltage().value(), frequency));
        preparedOperatingPoint.get_mutable_excitations_per_winding()[0] = excitation;
        make_waveform_size_power_of_two(&preparedOperatingPoint);
        auto preparedVoltage = preparedOperatingPoint.get_excitations_per_winding()[0].get_voltage().value();
        auto sampledVoltageWaveform = calculate_sampled_waveform(preparedVoltage.get_waveform().value(), frequency);

        std::string operatingPointLabel = "operating point " + std::to_string(operatingPointIndex);
        if (operatingPoint.get_name()) {
            operatingPointLabel += " ('" + operatingPoint.get_name().value() + "')";
        }
        std::string windingName = primaryWindingName ? primaryWindingName.value()
                                                     : excitation.get_name().value_or("unnamed");
        check_volt_second_balance(sampledVoltageWaveform,
            "Adviser inputs are invalid: " + operatingPointLabel + ", winding 0 ('" + windingName + "')");
    }
}

std::optional<SignalDescriptor> Inputs::get_single_winding_magnetizing_current(const OperatingPointExcitation& excitation,
                                                                              size_t numberWindings) {
    if (numberWindings != 1 || !excitation.get_current()) {
        return std::nullopt;
    }
    auto magnetizingCurrent = excitation.get_current().value();
    if (!magnetizingCurrent.get_waveform()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Excitation '" + excitation.get_name().value_or("unnamed") +
                                    "': the current has no waveform, so it cannot be taken as the magnetizing current");
    }
    // The current's processed data keeps the fields the caller gave (e.g. only offset and peak to
    // peak); the magnetizing current needs all of them (peak included), as
    // calculate_magnetizing_current would have produced.
    auto currentWaveform = magnetizingCurrent.get_waveform().value();
    Waveform sampledCurrentWaveform = is_waveform_sampled(currentWaveform) ? currentWaveform
                                                                            : calculate_sampled_waveform(currentWaveform, excitation.get_frequency());
    magnetizingCurrent.set_processed(calculate_processed_data(magnetizingCurrent, sampledCurrentWaveform, true));
    return magnetizingCurrent;
}

SignalDescriptor Inputs::calculate_magnetizing_current(const OperatingPoint& operatingPoint,
                                                       size_t windingIndex,
                                                       double magnetizingInductance,
                                                       const std::vector<double>& turnsRatios,
                                                       bool compress,
                                                       std::optional<Waveform> voltageSampledWaveform) {
    // The one decision for a winding's magnetizing current (ABT #1467): Faraday on the winding's
    // voltage, plus the DC the currents give. process_operating_point, the core adviser's
    // pre-processing and the autocomplete paths all come through here, so they cannot drift apart.
    const auto excitations = operatingPoint.get_excitations_per_winding();
    if (windingIndex >= excitations.size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "calculate_magnetizing_current: winding " + std::to_string(windingIndex) +
                                    " of an operating point with " + std::to_string(excitations.size()) + " windings");
    }
    if (magnetizingInductance <= 0) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "magnetizingInductance cannot be zero or negative");
    }
    auto excitation = excitations[windingIndex];
    // A single-winding inductor's winding current IS its magnetizing current (ABT #1585): the
    // voltage is integrated only when the current is absent or several windings share the flux.
    if (turnsRatios.empty()) {
        auto singleWindingMagnetizingCurrent = get_single_winding_magnetizing_current(excitation, excitations.size());
        if (singleWindingMagnetizingCurrent) {
            return singleWindingMagnetizingCurrent.value();
        }
    }
    if (!voltageSampledWaveform) {
        if (!excitation.get_voltage()) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "calculate_magnetizing_current: winding " + std::to_string(windingIndex) +
                                        " has no voltage to derive its magnetizing current from (Faraday)");
        }
        auto voltageExcitation = standardize_waveform(excitation.get_voltage().value(), excitation.get_frequency());
        auto waveform = voltageExcitation.get_waveform().value();
        if (!is_waveform_sampled(waveform)) {
            waveform = calculate_sampled_waveform(waveform, excitation.get_frequency());
        }
        voltageSampledWaveform = waveform;
    }

    double dcCurrent = 0;
    if (excitations.size() > 1) {
        // Several windings: the same classification as include_dc_offset_into_magnetizing_current,
        // on one conduction map shared with the anchor. No winding with a current: a voltage-only
        // operating point, no DC from the currents. Some with and some without: throws (ABT #1814).
        auto conduction = calculate_winding_conduction(operatingPoint);
        if (conduction && windings_alternate(conduction.value())) {
            dcCurrent = calculate_magnetizing_current_dc_anchor(excitation, conduction.value(), windingIndex, voltageSampledWaveform.value(), magnetizingInductance);
        }
    }
    return integrate_magnetizing_current(excitation, voltageSampledWaveform.value(), magnetizingInductance, compress, dcCurrent);
}

bool is_continuously_conducting_power(OperatingPointExcitation excitation) {
    if (!excitation.get_current()) {
        // If only voltage is defined, we are in transformer mode
        return true;
    }
    if (!excitation.get_voltage()) {
        // If only current is defined, we are in inductor mode
        return false;
    }
    if (!excitation.get_current()->get_waveform()) {
        throw std::invalid_argument("Missing current waveform");
    }
    if (!excitation.get_voltage()->get_waveform()) {
        throw std::invalid_argument("Missing voltage waveform");
    }
    std::vector<double> currentWaveform;
    std::vector<double> voltageWaveform;

    if (excitation.get_current()->get_waveform()->get_data().size() != excitation.get_voltage()->get_waveform()->get_data().size()) {
        // Determine a common number of points for both waveforms to ensure they have the same size after resampling
        size_t currentSize = excitation.get_current()->get_waveform()->get_data().size();
        size_t voltageSize = excitation.get_voltage()->get_waveform()->get_data().size();
        size_t numberPointsForSampling = settings.get_inputs_number_points_sampled_waveforms();
        size_t maxSize = std::max(currentSize, voltageSize);
        if (maxSize > numberPointsForSampling) {
            // Check if maxSize is a power of 2
            if (maxSize > 0 && ((maxSize & (maxSize - 1)) == 0)) {
                numberPointsForSampling = maxSize;
            }
            else {
                numberPointsForSampling = round_up_size_to_power_of_2(maxSize);
            }
        }
        auto currentSampledWaveform = Inputs::calculate_sampled_waveform(excitation.get_current()->get_waveform().value(), excitation.get_frequency(), numberPointsForSampling);
        auto voltageSampledWaveform = Inputs::calculate_sampled_waveform(excitation.get_voltage()->get_waveform().value(), excitation.get_frequency(), numberPointsForSampling);
        currentWaveform = currentSampledWaveform.get_data();
        voltageWaveform = voltageSampledWaveform.get_data();
    }
    else{
        currentWaveform = excitation.get_current()->get_waveform()->get_data();
        voltageWaveform = excitation.get_voltage()->get_waveform()->get_data();
    }
    std::vector<double> powerWaveform;
    for (size_t i = 0; i < currentWaveform.size(); ++i) {
        powerWaveform.push_back(currentWaveform[i] * voltageWaveform[i]);
    }

    double peakPower = *max_element(powerWaveform.begin(), powerWaveform.end());
    size_t numberPowerPointsUnderThreshold = 0;
    for (size_t i = 0; i < powerWaveform.size(); ++i) {
        if (fabs(powerWaveform[i]) < peakPower * 0.01) { // hardcoded
            numberPowerPointsUnderThreshold++;
        }
    }

    if (numberPowerPointsUnderThreshold > powerWaveform.size() * 0.1) { // hardcoded
        return false;
    }
    else {
        return true;
    }
}

SignalDescriptor Inputs::calculate_magnetizing_current(OperatingPointExcitation& excitation,
                                                                Waveform voltageSampledWaveform,
                                                                double magnetizingInductance,
                                                                bool compress,
                                                                double dcCurrent) {
    // The induced current asked for this voltage: an unbalanced one is refused with its own message (ABT #1855).
    check_volt_second_balance(voltageSampledWaveform, "Excitation '" + excitation.get_name().value_or("unnamed") + "'");
    return integrate_magnetizing_current(excitation, voltageSampledWaveform, magnetizingInductance, compress, dcCurrent);
}

SignalDescriptor Inputs::integrate_magnetizing_current(OperatingPointExcitation& excitation,
                                                       Waveform voltageSampledWaveform,
                                                       double magnetizingInductance,
                                                       bool compress,
                                                       double dcCurrent) {
    if (magnetizingInductance <= 0) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "magnetizingInductance cannot be zero or negative");
    }


    SignalDescriptor magnetizingCurrentExcitation;
    // Faraday, for every excitation (ABT #1467): v = N dPhi/dt, so the magnetizing current is the
    // voltage's integral over the magnetizing inductance, plus the DC anchor the caller derived
    // (dcCurrent). A FLYBACK_PRIMARY / FLYBACK_SECONDARY / UNIPOLAR_TRIANGULAR current used to be
    // rebuilt instead as a triangle from its own offset, peak to peak and duty cycle, rising over D
    // and falling over the whole of 1 - D: a DCM flyback lost its dead time (h1 of N_p i_m 17 % low
    // and 29 degrees off the windings' net MMF), and the load current's shape stood in for the flux.
    //
    // The integration constant is always removed: it depends only on where the numerical
    // integration starts within the cycle and is physically meaningless; any real DC bias is
    // supplied separately via dcCurrent. A small net volt-second imbalance (probe offset, resistive
    // drop, a capture that is not a whole number of periods) is removed before integrating.
    double meanVoltage = calculate_steady_state_mean_voltage(excitation, voltageSampledWaveform, "Excitation '" + excitation.get_name().value_or("unnamed") + "'");
    Waveform sampledMagnetizingCurrentWaveform = calculate_integral_waveform(sum_waveform(voltageSampledWaveform, -meanVoltage), true);
    sampledMagnetizingCurrentWaveform = multiply_waveform(sampledMagnetizingCurrentWaveform, 1.0 / magnetizingInductance);
    sampledMagnetizingCurrentWaveform = sum_waveform(sampledMagnetizingCurrentWaveform, dcCurrent);

    if (compress){
        // With the frequency: the sampled waveform's axis is closed (it reaches t0 + T), and closing
        // it again made the compressed magnetizing current span T + dt, which the sampler then
        // refuses (ABT #1837).
        magnetizingCurrentExcitation.set_waveform(compress_waveform(sampledMagnetizingCurrentWaveform, excitation.get_frequency()));
    }
    else {
        magnetizingCurrentExcitation.set_waveform(sampledMagnetizingCurrentWaveform);
    }

    magnetizingCurrentExcitation.set_harmonics(
        calculate_harmonics_data(sampledMagnetizingCurrentWaveform, sampledMagnetizingCurrentWaveform, excitation.get_frequency()));
    {
        auto processedData = calculate_processed_data(magnetizingCurrentExcitation, sampledMagnetizingCurrentWaveform);
        // The 'offset' field represents the AC mean of the magnetizing current (DC bias is captured in harmonics[0]).
        // We need to calculate the actual midpoint of the waveform (not rely on label-based offset which may return 0).
        // After subtracting dcCurrent, the result should be the offset of the underlying AC component.
        double waveformMidpoint = (*max_element(sampledMagnetizingCurrentWaveform.get_data().begin(), 
                                                 sampledMagnetizingCurrentWaveform.get_data().end()) +
                                   *min_element(sampledMagnetizingCurrentWaveform.get_data().begin(), 
                                                 sampledMagnetizingCurrentWaveform.get_data().end())) / 2;
        processedData.set_offset(waveformMidpoint - dcCurrent);
        magnetizingCurrentExcitation.set_processed(processedData);
    }

    return magnetizingCurrentExcitation;
}

SignalDescriptor Inputs::calculate_magnetizing_current(OperatingPointExcitation& excitation,
                                                                double magnetizingInductance,
                                                                bool compress,
                                                                double dcCurrent) {
    if (!excitation.get_voltage()) {
        throw std::invalid_argument("Missing voltage signal");
    }
    auto voltageExcitation = excitation.get_voltage().value();
    voltageExcitation = standardize_waveform(voltageExcitation, excitation.get_frequency());
    auto waveform = voltageExcitation.get_waveform().value();

    auto voltageSampledWaveform = calculate_sampled_waveform(waveform, excitation.get_frequency());
    return calculate_magnetizing_current(excitation, voltageSampledWaveform, magnetizingInductance, compress, dcCurrent);
}

SignalDescriptor Inputs::integrate_magnetizing_current(OperatingPointExcitation& excitation,
                                                       double magnetizingInductance,
                                                       bool compress,
                                                       double dcCurrent) {
    if (!excitation.get_voltage()) {
        throw std::invalid_argument("Missing voltage signal");
    }
    auto voltageExcitation = excitation.get_voltage().value();
    voltageExcitation = standardize_waveform(voltageExcitation, excitation.get_frequency());
    auto waveform = voltageExcitation.get_waveform().value();

    auto voltageSampledWaveform = calculate_sampled_waveform(waveform, excitation.get_frequency());
    return integrate_magnetizing_current(excitation, voltageSampledWaveform, magnetizingInductance, compress, dcCurrent);
}

OperatingPoint Inputs::process_operating_point(OperatingPoint operatingPoint, double magnetizingInductance, std::optional<std::vector<double>> turnsRatios, bool isDmcTopology,
                                                std::optional<std::vector<IsolationSide>> isolationSides) {
    // The stated frequency is the excitation's: it is never replaced by a harmonic of the waveform. A
    // line-cycle (PFC) record keeps its line frequency; its switching cycles are priced per cycle by
    // the core-loss models. A missing or invalid frequency throws (ABT #1467).
    for (size_t i = 0; i < operatingPoint.get_excitations_per_winding().size(); ++i) {
        double frequency = operatingPoint.get_excitations_per_winding()[i].get_frequency();
        if (!std::isfinite(frequency) || frequency <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "process_operating_point: excitation " + std::to_string(i) + " has frequency " + std::to_string(frequency) +
                "; every excitation needs its frequency (finite and positive)");
        }
    }

    std::vector<OperatingPointExcitation> processedExcitationsPerWinding;
    std::vector<Waveform> voltageSampledWaveforms;
    bool allExcitationHaveVoltageOrAlreadyCalculated = true;

    for (size_t windingIndex = 0; windingIndex < operatingPoint.get_excitations_per_winding().size(); ++windingIndex) {
        auto excitation = operatingPoint.get_excitations_per_winding()[windingIndex];
        // Here we processed this excitation current

        if (excitation.get_current()) {
            auto currentExcitation = excitation.get_current().value();

            if (!is_standardized(currentExcitation)) {
                currentExcitation = standardize_waveform(currentExcitation, excitation.get_frequency());
            }
            auto waveform = currentExcitation.get_waveform().value();
            Waveform sampledWaveform;
            if (!is_waveform_sampled(waveform)) {
                sampledWaveform = calculate_sampled_waveform(waveform, excitation.get_frequency());
            }
            else {
                sampledWaveform = waveform;
            }
            currentExcitation.set_harmonics(calculate_harmonics_data(waveform, sampledWaveform, excitation.get_frequency()));
            currentExcitation.set_processed(calculate_processed_data(currentExcitation, sampledWaveform, true, currentExcitation.get_processed()));
            excitation.set_current(currentExcitation);
        }
        else {
            if (operatingPoint.get_excitations_per_winding().size() == 2 && windingIndex == 1 && operatingPoint.get_excitations_per_winding()[0].get_current() && turnsRatios) {
                auto turnsRatio = turnsRatios.value()[0];
                auto currentExcitation = reflect_waveform(operatingPoint.get_excitations_per_winding()[0].get_current().value(), turnsRatio);
                auto waveform = currentExcitation.get_waveform().value();
                Waveform sampledWaveform;
                if (!is_waveform_sampled(waveform)) {
                    sampledWaveform = calculate_sampled_waveform(waveform, excitation.get_frequency());
                }
                else {
                    sampledWaveform = waveform;
                }
                currentExcitation.set_harmonics(calculate_harmonics_data(waveform, sampledWaveform, excitation.get_frequency()));
                currentExcitation.set_processed(calculate_processed_data(currentExcitation, sampledWaveform, true, currentExcitation.get_processed()));
                excitation.set_current(currentExcitation);
            }

        }
        // Here we processed this excitation voltage
        if (excitation.get_voltage()) {
            auto voltageExcitation = excitation.get_voltage().value();
            if (!is_standardized(voltageExcitation)) {
                voltageExcitation = standardize_waveform(voltageExcitation, excitation.get_frequency());
            }
            auto waveform = voltageExcitation.get_waveform().value();
            Waveform sampledWaveform;
            if (!is_waveform_sampled(waveform)) {
                sampledWaveform = calculate_sampled_waveform(waveform, excitation.get_frequency());
            }
            else {
                sampledWaveform = waveform;
            }
            voltageSampledWaveforms.push_back(sampledWaveform);
            voltageExcitation.set_harmonics(calculate_harmonics_data(waveform, sampledWaveform, excitation.get_frequency()));
            voltageExcitation.set_processed(calculate_processed_data(voltageExcitation, sampledWaveform));
            excitation.set_voltage(voltageExcitation);
            // A single-winding inductor's winding current IS its magnetizing current, exactly as in
            // the no-voltage branch below. Integrating the voltage (the generic path further down) is
            // needed only when the current is absent or several windings share the flux, and it was
            // wrong whenever the given voltage is not volt-second balanced: Kirchhoff's line-cycle
            // PFC voltages integrated to kA "magnetizing" peaks against ~12 A of real current, and
            // no core passed the energy-stored gate (ABT #1585).
            if (!excitation.get_magnetizing_current()) {
                auto singleWindingMagnetizingCurrent = get_single_winding_magnetizing_current(excitation, operatingPoint.get_excitations_per_winding().size());
                if (singleWindingMagnetizingCurrent) {
                    excitation.set_magnetizing_current(singleWindingMagnetizingCurrent.value());
                }
            }
        }
        else {
            if (operatingPoint.get_excitations_per_winding().size() == 1 && excitation.get_current()) {
                auto voltageExcitation = calculate_induced_voltage(excitation, magnetizingInductance);
                excitation.set_voltage(voltageExcitation);
                excitation.set_magnetizing_current(excitation.get_current().value());
                allExcitationHaveVoltageOrAlreadyCalculated = false;
            }
            else if (operatingPoint.get_excitations_per_winding().size() == 2 && windingIndex == 1 && operatingPoint.get_excitations_per_winding()[0].get_voltage() && turnsRatios) {
                auto turnsRatio = turnsRatios.value()[0];
                auto voltageExcitation = reflect_waveform(operatingPoint.get_excitations_per_winding()[0].get_voltage().value(), 1 / turnsRatio);
                auto waveform = voltageExcitation.get_waveform().value();
                Waveform sampledWaveform;
                if (!is_waveform_sampled(waveform)) {
                    sampledWaveform = calculate_sampled_waveform(waveform, excitation.get_frequency());
                }
                else {
                    sampledWaveform = waveform;
                }
                voltageSampledWaveforms.push_back(sampledWaveform);
                voltageExcitation.set_harmonics(calculate_harmonics_data(waveform, sampledWaveform, excitation.get_frequency()));
                voltageExcitation.set_processed(calculate_processed_data(voltageExcitation, sampledWaveform, true, voltageExcitation.get_processed()));
                excitation.set_voltage(voltageExcitation);
            }
            else {
                allExcitationHaveVoltageOrAlreadyCalculated = false;
            }
        }
        processedExcitationsPerWinding.push_back(excitation);
    }
    operatingPoint.set_excitations_per_winding(processedExcitationsPerWinding);

    if (allExcitationHaveVoltageOrAlreadyCalculated) {

        std::vector<double> turnsRatios;
        double primaryVoltageRms = operatingPoint.get_excitations_per_winding()[0].get_voltage()->get_processed()->get_rms().value();
        for (size_t windingIndex = 1; windingIndex < operatingPoint.get_excitations_per_winding().size(); ++windingIndex) {
            auto excitation = operatingPoint.get_excitations_per_winding()[windingIndex];
            double turnsRatio = primaryVoltageRms / excitation.get_voltage()->get_processed()->get_rms().value();

            turnsRatios.push_back(turnsRatio);
        }

        // ── CMC handling ─────────────────────────────────────────────────
        // In a common-mode choke every winding carries identical DM + CM
        // currents (MKF stores every winding's current "into the dot"), so
        // the DM components cancel in the net core MMF and only the CM
        // ripple drives flux. The generic per-winding path below integrates
        // each winding's voltage and adds that winding's DC offset back in,
        // which for a CMC means N · (I_dm_dc + I_cm_ac) instead of just
        // N · I_cm_ac — over-counts B by a factor of (I_line / I_cm) and
        // can make the core look saturated when it isn't.
        //
        // When can_be_common_mode_choke(...) fires (2+ matched winding
        // currents at identical frequencies and amplitudes) plug in the
        // CM-only magnetizing current across every winding — same logic
        // MagnetizingInductance.cpp:149 already uses before core losses run.
        //
        // Guard: a 1:1 transformer test scenario looks the same to
        // can_be_common_mode_choke (symmetric currents and voltages), so
        // only fire when the caller has NOT already pre-populated a
        // magnetizing_current. Real CMC paths leave it empty here; test
        // scenarios (create_quick_operating_point with turnsRatios={1})
        // pre-populate it to pin a specific answer, and stomping that was
        // regressing TestCoreAdviserAvailableCores on 1:1 transformers.
        bool primaryHasMagnetizingCurrent =
            !operatingPoint.get_excitations_per_winding().empty()
            && operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current().has_value();
        // ── DMC handling ─────────────────────────────────────────────────
        // Multi-winding DMC topologies (SINGLE_PHASE_BALANCED, THREE_PHASE,
        // THREE_PHASE_WITH_NEUTRAL) carry one phase-current per winding in
        // physical sign convention. The flux-driving "magnetizing" current
        // the core sees is a case-specific combination (see
        // get_differential_mode_choke_magnetizing_current). Stamp the
        // synthesized current onto winding[0] so MagnetizingInductance and
        // the downstream loss/flux pipeline pick it up via the standard
        // pre-stamped path.
        // Single-phase 1-winding DMC is intentionally NOT intercepted here
        // — it's a plain series inductor and the regular path handles it.
        bool isDmcOperatingPoint = !primaryHasMagnetizingCurrent
            && can_be_differential_mode_choke(operatingPoint, isDmcTopology);
        // The windings' isolation sides decide whether identical currents add (one side: a
        // common-mode choke) or cancel (a transformer). Undeclared, MKF's index convention applies.
        std::vector<IsolationSide> windingIsolationSides;
        // (Fewer declared than windings: the rest by index, as get_isolation_sides_used extends them.)
        if (isolationSides) {
            windingIsolationSides = isolationSides.value();
        }
        for (size_t windingIndex = windingIsolationSides.size(); windingIndex < operatingPoint.get_excitations_per_winding().size(); ++windingIndex) {
            windingIsolationSides.push_back(get_isolation_side_from_index(windingIndex));
        }
        bool isCmcOperatingPoint = !primaryHasMagnetizingCurrent
            && !isDmcOperatingPoint
            && can_be_common_mode_choke(operatingPoint, windingIsolationSides);
        if (isDmcOperatingPoint) {
            auto dmMagnetizingCurrent = get_differential_mode_choke_magnetizing_current(operatingPoint);
            // Stamp the synthesized DMC magnetizing current on every winding;
            // downstream consumers read winding[0] but mirroring keeps the
            // operating point internally consistent (each winding "sees" the
            // same coupled flux).
            for (size_t windingIndex = 0; windingIndex < processedExcitationsPerWinding.size(); ++windingIndex) {
                processedExcitationsPerWinding[windingIndex].set_magnetizing_current(dmMagnetizingCurrent);
            }
        }
        else if (isCmcOperatingPoint) {
            auto cmMagnetizingCurrent = get_common_mode_choke_magnetizing_current(operatingPoint);
            for (size_t windingIndex = 0; windingIndex < processedExcitationsPerWinding.size(); ++windingIndex) {
                processedExcitationsPerWinding[windingIndex].set_magnetizing_current(cmMagnetizingCurrent);
            }
        }
        else {
            for (size_t windingIndex = 0; windingIndex < operatingPoint.get_excitations_per_winding().size(); ++windingIndex) {
                auto excitation = operatingPoint.get_excitations_per_winding()[windingIndex];

                if (!excitation.get_magnetizing_current() && magnetizingInductance > 0) {
                    Waveform waveform = excitation.get_voltage()->get_waveform().value();
                    if (windingIndex < voltageSampledWaveforms.size()) {
                        waveform = voltageSampledWaveforms[windingIndex];
                    }
                    excitation.set_magnetizing_current(calculate_magnetizing_current(operatingPoint, windingIndex, magnetizingInductance, turnsRatios, false, waveform));
                }
                processedExcitationsPerWinding[windingIndex] = excitation;
            }
        }

    }
    operatingPoint.set_excitations_per_winding(processedExcitationsPerWinding);
    return operatingPoint;
}

void Inputs::process(std::optional<std::variant<double, std::vector<double>>> magnetizingInductance) {
    auto operatingPoints = get_mutable_operating_points();
    std::vector<OperatingPoint> processed_operating_points;

    // Detect DMC topology from design requirements so process_operating_point
    // can synthesize the effective core magnetizing current for multi-winding
    // DMCs (configuration is implied by the operating point's winding count).
    bool isDmcTopology = get_design_requirements().get_topology().has_value()
        && get_design_requirements().get_topology().value() == MAS::Topology::DIFFERENTIAL_MODE_CHOKE;

    // The windings' isolation sides (ABT #1467): declared, or every winding on the line side for a
    // declared common-mode choke; otherwise process_operating_point applies MKF's index convention.
    std::optional<std::vector<IsolationSide>> isolationSides = get_design_requirements().get_isolation_sides();
    bool isCmcTopology = get_design_requirements().get_topology().has_value()
        && get_design_requirements().get_topology().value() == MAS::Topology::COMMON_MODE_CHOKE;
    if (!isolationSides && isCmcTopology) {
        size_t numberWindings = get_design_requirements().get_turns_ratios().size() + 1;
        for (auto& operatingPoint : operatingPoints) {
            numberWindings = std::max(numberWindings, operatingPoint.get_excitations_per_winding().size());
        }
        isolationSides = std::vector<IsolationSide>(numberWindings, IsolationSide::PRIMARY);
    }

    for (size_t operatingPointIndex = 0; operatingPointIndex < operatingPoints.size(); ++operatingPointIndex) {
        auto operatingPoint = operatingPoints[operatingPointIndex];
        double magnetizingInductanceToProcess;
        if (magnetizingInductance) {
            if (std::holds_alternative<std::vector<double>>(magnetizingInductance.value())) {
                auto magnetizingInductancePerPoint = std::get<std::vector<double>>(magnetizingInductance.value());
                if (magnetizingInductancePerPoint.size() == 0) {
                    magnetizingInductanceToProcess = resolve_dimensional_values(get_design_requirements().get_magnetizing_inductance());
                }
                else if (operatingPointIndex < magnetizingInductancePerPoint.size()) {
                    magnetizingInductanceToProcess = magnetizingInductancePerPoint[operatingPointIndex];
                }
                else {
                    magnetizingInductanceToProcess = magnetizingInductancePerPoint.back();
                }
            }
            else {
                magnetizingInductanceToProcess = std::get<double>(magnetizingInductance.value());
            }
        }
        else {
            magnetizingInductanceToProcess = resolve_dimensional_values(get_design_requirements().get_magnetizing_inductance());
        }
        processed_operating_points.push_back(process_operating_point(operatingPoint, magnetizingInductanceToProcess, std::nullopt, isDmcTopology, isolationSides));
    }
    set_operating_points(processed_operating_points);
}


OperatingPoint Inputs::create_operating_point_with_sinusoidal_current_mask(double frequency,
                                                                                  double magnetizingInductance,
                                                                                  double temperature,
                                                                                  std::vector<double> turnsRatios,
                                                                                  std::vector<double> currentPeakMask,
                                                                                  double currentOffset) {

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
        ProcessedWaveform processed;
        processed.set_label(WaveformLabel::SINUSOIDAL);
        processed.set_peak_to_peak(2 * currentPeakMask[0]);
        processed.set_duty_cycle(0.5);
        processed.set_offset(currentOffset);
        current.set_processed(processed);
        current = standardize_waveform(current, frequency);
        current.set_processed(calculate_basic_processed_data(current.get_waveform().value()));
        excitation.set_current(current);
        if (magnetizingInductance > 0) {
            auto voltage = calculate_induced_voltage(excitation, magnetizingInductance);
            voltage.set_processed(calculate_basic_processed_data(voltage.get_waveform().value()));
            excitation.set_voltage(voltage);
            auto magnetizingCurrent = integrate_magnetizing_current(excitation, magnetizingInductance, true, 0.0);
            magnetizingCurrent.set_processed(calculate_basic_processed_data(magnetizingCurrent.get_waveform().value()));
            excitation.set_magnetizing_current(magnetizingCurrent);
        }
        operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
    }
    for (size_t turnsRatioIndex = 0; turnsRatioIndex < turnsRatios.size(); ++turnsRatioIndex) {
        auto turnsRatio = turnsRatios[turnsRatioIndex];
        OperatingPointExcitation excitation;
        excitation.set_frequency(frequency);
        SignalDescriptor current;
        ProcessedWaveform processed;
        processed.set_label(WaveformLabel::SINUSOIDAL);
        processed.set_peak_to_peak(2 * currentPeakMask[turnsRatioIndex + 1]);
        processed.set_duty_cycle(0.5);
        processed.set_offset(0);
        current.set_processed(processed);
        current = standardize_waveform(current, frequency);
        excitation.set_current(current);
        if (magnetizingInductance > 0) {
            auto voltage = calculate_induced_voltage(excitation, magnetizingInductance / pow(turnsRatio, 2));
            voltage.set_processed(calculate_basic_processed_data(voltage.get_waveform().value()));
            excitation.set_voltage(voltage);
            auto magnetizingCurrent = integrate_magnetizing_current(excitation, magnetizingInductance, true, 0.0);
            magnetizingCurrent.set_processed(calculate_basic_processed_data(magnetizingCurrent.get_waveform().value()));
            excitation.set_magnetizing_current(magnetizingCurrent);
        }
        operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
    }

    return process_operating_point(operatingPoint, magnetizingInductance);

}

Inputs Inputs::create_quick_operating_point(double frequency,
                                                          double magnetizingInductance,
                                                          double temperature,
                                                          WaveformLabel waveShape,
                                                          double peakToPeak,
                                                          double dutyCycle,
                                                          double dcCurrent,
                                                          std::vector<double> turnsRatios) {

    Inputs inputs;

    DesignRequirements designRequirements;
    DimensionWithTolerance magnetizingInductanceWithTolerance;
    magnetizingInductanceWithTolerance.set_minimum(magnetizingInductance * 0.8);
    magnetizingInductanceWithTolerance.set_nominal(magnetizingInductance);
    magnetizingInductanceWithTolerance.set_maximum(magnetizingInductance * 1.2);

    InsulationRequirements insulationRequirements;
    auto overvoltageCategory = OvervoltageCategory::II;
    auto cti = Cti::GROUP_I;
    DimensionWithTolerance altitude;
    DimensionWithTolerance mainSupplyVoltage;
    auto pollutionDegree = PollutionDegree::PD1;
    auto standards = std::vector<InsulationStandards>{};
    altitude.set_maximum(2000);
    mainSupplyVoltage.set_nominal(400);
    auto insulationType = IsolationClass::BASIC;

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
        ProcessedWaveform processed;
        processed.set_label(waveShape);
        processed.set_peak_to_peak(peakToPeak);
        processed.set_duty_cycle(dutyCycle);
        processed.set_offset(0);
        voltage.set_processed(processed);
        voltage = standardize_waveform(voltage, frequency);

        excitation.set_voltage(voltage);
        if (magnetizingInductance > 0) {
            auto current = integrate_magnetizing_current(excitation, magnetizingInductance, true, dcCurrent);
            excitation.set_current(current);
            excitation.set_magnetizing_current(current);
        }
        operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
    }
    for (auto& turnsRatio : turnsRatios) {
        OperatingPointExcitation excitation;
        excitation.set_frequency(frequency);
        SignalDescriptor voltage;
        ProcessedWaveform processed;
        processed.set_label(waveShape);
        processed.set_peak_to_peak(peakToPeak * turnsRatio);
        processed.set_duty_cycle(dutyCycle);
        processed.set_offset(0);
        voltage.set_processed(processed);
        voltage = standardize_waveform(voltage, frequency);
        excitation.set_voltage(voltage);
        if (magnetizingInductance > 0) {
            auto current = integrate_magnetizing_current(excitation, magnetizingInductance * pow(turnsRatio, 2), true, dcCurrent);
            excitation.set_current(current);
            excitation.set_magnetizing_current(current);
        }
        operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
    }

    inputs.get_mutable_operating_points().push_back(operatingPoint);
    inputs.process();
    return inputs;

}
Inputs Inputs::create_quick_operating_point_only_current(double frequency,
                                                          double magnetizingInductance,
                                                          double temperature,
                                                          WaveformLabel waveShape,
                                                          double peakToPeak,
                                                          double dutyCycle,
                                                          double dcCurrent,
                                                          std::vector<double> turnsRatios) {

    Inputs inputs;

    DesignRequirements designRequirements;
    DimensionWithTolerance magnetizingInductanceWithTolerance;
    magnetizingInductanceWithTolerance.set_minimum(magnetizingInductance * 0.8);
    magnetizingInductanceWithTolerance.set_nominal(magnetizingInductance);
    magnetizingInductanceWithTolerance.set_maximum(magnetizingInductance * 1.2);
    InsulationRequirements insulationRequirements;

    auto overvoltageCategory = OvervoltageCategory::II;
    auto cti = Cti::GROUP_I;
    DimensionWithTolerance altitude;
    DimensionWithTolerance mainSupplyVoltage;
    auto pollutionDegree = PollutionDegree::PD1;
    auto standards = std::vector<InsulationStandards>{};
    altitude.set_maximum(2000);
    mainSupplyVoltage.set_nominal(400);
    auto insulationType = IsolationClass::BASIC;

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
        ProcessedWaveform processed;
        processed.set_label(waveShape);
        processed.set_peak(peakToPeak / 2 + dcCurrent);
        processed.set_peak_to_peak(peakToPeak);
        processed.set_duty_cycle(dutyCycle);
        processed.set_offset(dcCurrent);
        current.set_processed(processed);
        current = standardize_waveform(current, frequency);
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
        ProcessedWaveform processed;
        processed.set_label(waveShape);
        processed.set_peak(peakToPeak / 2 * turnsRatio);
        processed.set_peak_to_peak(peakToPeak * turnsRatio);
        processed.set_duty_cycle(dutyCycle);
        processed.set_offset(dcCurrent);
        current.set_processed(processed);
        current = standardize_waveform(current, frequency);
        excitation.set_current(current);
        if (magnetizingInductance > 0) {
            auto voltage = calculate_induced_voltage(excitation, magnetizingInductance / pow(turnsRatio, 2));
            excitation.set_voltage(voltage);
        }
        operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
    }

    inputs.get_mutable_operating_points().push_back(operatingPoint);
    inputs.process();

    return inputs;
}
Inputs Inputs::create_quick_operating_point_only_current(double frequency,
                                                          double magnetizingInductance,
                                                          double temperature,
                                                          WaveformLabel waveShape,
                                                          std::vector<double> peakToPeaks,
                                                          double dutyCycle,
                                                          double dcCurrent,
                                                          std::vector<double> turnsRatios) {

    Inputs inputs;

    DesignRequirements designRequirements;
    DimensionWithTolerance magnetizingInductanceWithTolerance;
    magnetizingInductanceWithTolerance.set_minimum(magnetizingInductance * 0.8);
    magnetizingInductanceWithTolerance.set_nominal(magnetizingInductance);
    magnetizingInductanceWithTolerance.set_maximum(magnetizingInductance * 1.2);
    InsulationRequirements insulationRequirements;

    auto overvoltageCategory = OvervoltageCategory::II;
    auto cti = Cti::GROUP_I;
    DimensionWithTolerance altitude;
    DimensionWithTolerance mainSupplyVoltage;
    auto pollutionDegree = PollutionDegree::PD1;
    auto standards = std::vector<InsulationStandards>{};
    altitude.set_maximum(2000);
    mainSupplyVoltage.set_nominal(400);
    auto insulationType = IsolationClass::BASIC;

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
        ProcessedWaveform processed;
        processed.set_label(waveShape);
        processed.set_peak(peakToPeak / 2 + dcCurrent);
        processed.set_peak_to_peak(peakToPeak);
        processed.set_duty_cycle(dutyCycle);
        processed.set_offset(dcCurrent);
        current.set_processed(processed);
        current = standardize_waveform(current, frequency);
        excitation.set_current(current);
        if (magnetizingInductance > 0) {
            auto voltage = calculate_induced_voltage(excitation, magnetizingInductance / pow(turnsRatio, 2));
            excitation.set_voltage(voltage);
        }
        operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
    }

    inputs.get_mutable_operating_points().push_back(operatingPoint);
    inputs.process();

    return inputs;
}

OperatingPoint Inputs::get_operating_point(size_t index) {
    return get_mutable_operating_points()[index];
}

OperatingPointExcitation Inputs::get_winding_excitation(size_t operatingPointIndex, size_t windingIndex) {
    return get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex];
}

OperatingPointExcitation Inputs::get_primary_excitation(size_t operatingPointIndex) {
    return get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[0];
}

OperatingPointExcitation Inputs::get_primary_excitation(const OperatingPoint& operatingPoint) {
    return operatingPoint.get_excitations_per_winding()[0];
}

void Inputs::make_waveform_size_power_of_two(OperatingPoint* operatingPoint) {
    double frequency = operatingPoint->get_excitations_per_winding()[0].get_frequency();

    // Phase 5 fix: iterate over ALL windings, not just the primary.
    // CMCs (and any multi-winding design where can_be_common_mode_choke
    // matters) require every winding's current to carry identical harmonics
    // after standardization.  Only touching winding 0 left the secondary
    // with its original (possibly different) harmonics, causing
    // can_be_common_mode_choke to return false and routing the CMC
    // through the single-winding path that uses the full line current as
    // magnetizing current — producing false saturation.
    for (size_t w = 0; w < operatingPoint->get_mutable_excitations_per_winding().size(); ++w) {
        OperatingPointExcitation excitation = operatingPoint->get_excitations_per_winding()[w];

        // Skip windings whose current/voltage carries only `processed`
        // (RMS/peak) without an explicit waveform — there is nothing to
        // resample, and `make_waveform_size_power_of_two` has no business
        // failing the caller in that case. Anything downstream that *does*
        // need a populated waveform will throw at its own call site with a
        // meaningful message rather than a bare `bad_optional_access` here.
        if (excitation.get_current()) {
            auto current = operatingPoint->get_mutable_excitations_per_winding()[w].get_current().value();
            if (current.get_waveform()) {
                auto currentWaveform = current.get_waveform().value();
                if (!is_size_power_of_2(currentWaveform.get_data())) {
                    auto currentSampledWaveform = Inputs::calculate_sampled_waveform(currentWaveform, frequency);
                    current.set_waveform(currentSampledWaveform);
                    current.set_harmonics(calculate_harmonics_data(currentWaveform, currentSampledWaveform, frequency));
                    current.set_processed(calculate_processed_data(current, currentSampledWaveform, true, current.get_processed()));
                    operatingPoint->get_mutable_excitations_per_winding()[w].set_current(current);
                }
            }
        }
        if (excitation.get_voltage()) {
            auto voltage = operatingPoint->get_mutable_excitations_per_winding()[w].get_voltage().value();
            if (voltage.get_waveform()) {
                auto voltageWaveform = voltage.get_waveform().value();
                if (!is_size_power_of_2(voltageWaveform.get_data())) {
                    auto voltageSampledWaveform = Inputs::calculate_sampled_waveform(voltageWaveform, frequency);
                    voltage.set_waveform(voltageSampledWaveform);
                    operatingPoint->get_mutable_excitations_per_winding()[w].set_voltage(voltage);
                }
            }
        }
    }
}

double Inputs::calculate_waveform_coefficient(OperatingPoint* operatingPoint) {
    OperatingPointExcitation excitation = Inputs::get_primary_excitation(*operatingPoint);
    double frequency = excitation.get_frequency();
    Waveform sampledWaveform = excitation.get_voltage()->get_waveform().value();

    if (sampledWaveform.get_time() && sampledWaveform.get_data().size() < settings.get_inputs_number_points_sampled_waveforms()) {
        sampledWaveform = calculate_sampled_waveform(sampledWaveform, frequency);
    }

    std::vector<double> source = sampledWaveform.get_data();
    std::vector<double> integration;

    double integral = 0;
    integration.push_back(integral);
    double time_per_point = 1 / frequency / source.size();
    source = std::vector<double>(source.begin(), source.end() - source.size() / 2);

    for (size_t i = 0; i < source.size() - 1; i++) {
        // Trapezoidal rule: (s[i] + s[i+1]) / 2. The previous fabs() on the
        // difference inflated the integral on every falling segment.
        integral += (source[i + 1] + source[i]) / 2.;
    }
    integral *= time_per_point;

    double voltageRms = excitation.get_voltage()->get_processed()->get_rms().value();

    double waveformCoefficient = 2 * voltageRms / (frequency * integral);
    return waveformCoefficient;
}

double Inputs::calculate_instantaneous_power(OperatingPointExcitation excitation) {
    double frequency = excitation.get_frequency();
    if (!excitation.get_voltage()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Voltage signal is missing");
    }
    if (!excitation.get_voltage()->get_waveform()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Voltage waveform is missing");
    }
    if (!excitation.get_current()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Current signal is missing");
    }
    if (!excitation.get_current()->get_waveform()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Current waveform is missing");
    }
    Waveform voltageSampledWaveform = excitation.get_voltage()->get_waveform().value();
    Waveform currentSampledWaveform = excitation.get_current()->get_waveform().value();

    // Both signals go on one grid, as fine as the finer of them, and |v*i| is averaged over ALL of
    // it. The loop used to stop at the 128-sample setting while the sampler keeps a denser
    // waveform's own resolution, so an imported waveform (8192 samples) was averaged over the
    // first 1.6 % of its period and the figure moved with the file's resolution (ABT #1325).
    size_t numberPoints = settings.get_inputs_number_points_sampled_waveforms();
    for (const auto* waveform : {&voltageSampledWaveform, &currentSampledWaveform}) {
        size_t size = waveform->get_data().size();
        if (size > numberPoints) {
            // Exact powers of 2 are kept as they are: the rounding goes through a floating log.
            numberPoints = (size & (size - 1)) == 0? size : round_up_size_to_power_of_2(size);
        }
    }
    if (voltageSampledWaveform.get_data().size() != numberPoints) {
        voltageSampledWaveform = calculate_sampled_waveform(voltageSampledWaveform, frequency, numberPoints);
    }
    if (currentSampledWaveform.get_data().size() != numberPoints) {
        currentSampledWaveform = calculate_sampled_waveform(currentSampledWaveform, frequency, numberPoints);
    }
    if (voltageSampledWaveform.get_data().size() != currentSampledWaveform.get_data().size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "calculate_instantaneous_power: voltage and current sampled to " +
            std::to_string(voltageSampledWaveform.get_data().size()) + " and " +
            std::to_string(currentSampledWaveform.get_data().size()) + " points; they must share one grid");
    }

    std::vector<double> powerPoints;

    for (size_t i = 0; i < voltageSampledWaveform.get_data().size(); i++) {
        powerPoints.push_back(fabs(voltageSampledWaveform.get_data()[i] * currentSampledWaveform.get_data()[i]));
    }

    double instantaneousPower = std::accumulate(std::begin(powerPoints), std::end(powerPoints), 0.0) /
                             powerPoints.size();
    return instantaneousPower;
}

WaveformLabel Inputs::try_guess_waveform_label(Waveform waveform) {
    // Delegates to WaveformProcessor, which fixed this classifier (ABT #602):
    // the sine-vs-custom fallback divided its error metric by N twice
    // (error /= data.size() then error /= area, where area was already a sum
    // over N points), so the 5% acceptance threshold became an effective
    // 2560% relative error at N=512 and every imported waveform — rectangles
    // and triangles included — was labelled sinusoidal; the reference sine
    // also ran off the nominal sample count instead of this waveform's own.
    // This used to be a verbatim, independently-bugged twin of that logic —
    // kept as a delegating wrapper so the two can no longer drift apart,
    // matching every other Inputs::X waveform-DSP method in this file.
    return WaveformProcessor::try_guess_waveform_label(waveform, settings.get_inputs_number_points_sampled_waveforms());
}

void Inputs::scale_time_to_frequency(Inputs& inputs, double newFrequency, bool cleanFrequencyDependentFields, bool processSignals, bool useCurrentAsBase){
    for (auto& operatingPoint : inputs.get_mutable_operating_points()) {
        Inputs::scale_time_to_frequency(operatingPoint, newFrequency, cleanFrequencyDependentFields, processSignals, useCurrentAsBase);
    }}

void Inputs::scale_time_to_frequency(OperatingPoint& operatingPoint, double newFrequency, bool cleanFrequencyDependentFields, bool processSignals, bool useCurrentAsBase){
    for (auto& excitation : operatingPoint.get_mutable_excitations_per_winding()) {
        scale_time_to_frequency(excitation, newFrequency, cleanFrequencyDependentFields, processSignals, useCurrentAsBase);
    }}

void Inputs::scale_time_to_frequency(OperatingPointExcitation& excitation, double newFrequency, bool cleanFrequencyDependentFields, bool processSignals, bool useCurrentAsBase){
    excitation.set_frequency(newFrequency);
    // A signal whose time axis is rescaled keeps its waveform; its harmonics and processed data were
    // derived from the old time axis. They are recomputed from the rescaled waveform when processSignals
    // is set, and dropped otherwise, so that no consumer reads a spectrum at the old frequencies (a
    // swept DMC record kept its 100 kHz ripple harmonic while its waveform moved to 1 MHz, ABT #1764).
    // A signal without a waveform is described by its processed data alone and is left as it is.
    auto rescale = [&](SignalDescriptor signal) {
        signal.set_waveform(scale_time_to_frequency(signal.get_waveform().value(), newFrequency));
        if (processSignals) {
            auto sampledWaveform = Inputs::calculate_sampled_waveform(signal.get_waveform().value(), newFrequency);
            signal.set_harmonics(Inputs::calculate_harmonics_data(signal.get_waveform().value(), sampledWaveform, newFrequency));
            signal.set_processed(Inputs::calculate_processed_data(signal, sampledWaveform, true));
        }
        else {
            signal.set_harmonics(std::nullopt);
            signal.set_processed(std::nullopt);
        }
        return signal;
    };
    if (excitation.get_current() && excitation.get_current()->get_waveform()) {
        excitation.set_current(rescale(excitation.get_current().value()));
    }
    if (excitation.get_voltage() && excitation.get_voltage()->get_waveform()) {
        excitation.set_voltage(rescale(excitation.get_voltage().value()));
    }
    if (cleanFrequencyDependentFields) {
        if (useCurrentAsBase) {
            excitation.set_voltage(std::nullopt);
        }
        else {
            excitation.set_current(std::nullopt);
        }
        excitation.set_magnetizing_current(std::nullopt);
        excitation.set_magnetic_flux_density(std::nullopt);
        excitation.set_magnetic_field_strength(std::nullopt);
    }
    else {
        if (excitation.get_magnetizing_current() && excitation.get_magnetizing_current()->get_waveform()) {
            excitation.set_magnetizing_current(rescale(excitation.get_magnetizing_current().value()));
        }
        if (excitation.get_magnetic_flux_density() && excitation.get_magnetic_flux_density()->get_waveform()) {
            excitation.set_magnetic_flux_density(rescale(excitation.get_magnetic_flux_density().value()));
        }
        if (excitation.get_magnetic_field_strength() && excitation.get_magnetic_field_strength()->get_waveform()) {
            excitation.set_magnetic_field_strength(rescale(excitation.get_magnetic_field_strength().value()));
        }
    }
}
Waveform Inputs::scale_time_to_frequency(Waveform waveform, double newFrequency){
    std::vector<double> scaledTime;
    std::vector<double> time = waveform.get_time().value();
    double oldFrequency = 1.0 / (time.back() - time.front());
    for (auto& timePoint : time) {
        scaledTime.push_back(timePoint * oldFrequency / newFrequency);
    }
    waveform.set_time(scaledTime);
    return waveform;
}

void process_voltage(OperatingPointExcitation& excitation) {
    // ABT #1330: a voltage described by its processed parameters alone (label, peakToPeak, offset,
    // dutyCycle) is rebuilt with MKF's own create_waveform at the excitation's frequency, the same
    // reconstruction autocomplete's inputs step does, instead of refusing it. Nothing is assumed: a
    // processed block that cannot define the waveform throws naming what is missing (create_waveform
    // itself would default a missing duty cycle to 0.5 and a dead time to 0).
    auto voltage = excitation.get_voltage().value();
    if (!voltage.get_waveform()) {
        if (!voltage.get_processed()) {
            throw InvalidInputException(ErrorCode::MISSING_DATA,
                "The voltage has neither a waveform nor processed data to build one from");
        }
        Inputs::throw_if_processed_cannot_define_waveform(voltage.get_processed().value(), "voltage");
        voltage = Inputs::standardize_waveform(voltage, excitation.get_frequency());
    }
    ProcessedWaveform processed = Inputs::calculate_processed_data(voltage.get_waveform().value(), excitation.get_frequency());
    voltage.set_processed(processed);
    excitation.set_voltage(voltage);
}


double Inputs::get_maximum_voltage_peak() {
    if (get_operating_points().size() == 0)
        throw std::invalid_argument("There are no operating points");

    double maximumVoltage = 0;
    for (auto& operatingPoint : get_operating_points()) {
        if (operatingPoint.get_excitations_per_winding().size() == 0)
            throw std::invalid_argument("There are no winding excitation in operating point");

        for (auto excitation : operatingPoint.get_excitations_per_winding()) {
            if (!excitation.get_voltage()) {
                throw std::invalid_argument("Missing voltage in excitation");
            }
            if (!excitation.get_voltage()->get_processed()) {
                process_voltage(excitation);
                // throw std::invalid_argument("Voltage has not been processed");
            }
            else if (!excitation.get_voltage()->get_processed()->get_peak()) {
                process_voltage(excitation);
            }
            maximumVoltage = std::max(maximumVoltage, excitation.get_voltage()->get_processed()->get_peak().value());
        }
    }

    return maximumVoltage;
}

double Inputs::get_maximum_voltage_rms() {
    if (get_operating_points().size() == 0)
        throw std::invalid_argument("There are no operating points");

    double maximumVoltage = 0;
    for (auto& operatingPoint : get_mutable_operating_points()) {
        if (operatingPoint.get_excitations_per_winding().size() == 0)
            throw std::invalid_argument("There are no winding excitation in operating point");

        for (auto excitation : operatingPoint.get_mutable_excitations_per_winding()) {
            if (!excitation.get_voltage()) 
                throw std::invalid_argument("Missing voltage in excitation");
            if (!excitation.get_voltage()->get_processed()) {
                process_voltage(excitation);
                // throw std::invalid_argument("Voltage has not been processed");
            }
            if (!excitation.get_voltage()->get_processed()->get_rms()) {
                process_voltage(excitation);
                // ProcessedWaveform processed = calculate_processed_data(excitation.get_voltage()->get_waveform().value(), excitation.get_frequency());
                // auto voltage = excitation.get_voltage().value();
                // voltage.set_processed(processed);
                // excitation.set_voltage(voltage);
                // throw std::invalid_argument("RMS has not been calculated");
            }
            maximumVoltage = std::max(maximumVoltage, excitation.get_voltage()->get_processed()->get_rms().value());
        }
    }

    return maximumVoltage;
}

double Inputs::get_maximum_current_rms() {
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

double Inputs::get_maximum_current_peak() {
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

double Inputs::get_maximum_voltage_peak(size_t windingIndex) {
    if (get_operating_points().size() == 0)
        throw std::invalid_argument("There are no operating points");

    double maximumVoltage = 0;
    for (auto& operatingPoint : get_operating_points()) {
        if (operatingPoint.get_excitations_per_winding().size() == 0)
            throw std::invalid_argument("There are no winding excitation in operating point");

        auto excitation = operatingPoint.get_excitations_per_winding()[windingIndex];
        if (!excitation.get_voltage()) 
            throw std::invalid_argument("Missing voltage in excitation");
        if (!excitation.get_voltage()->get_processed()) {
            process_voltage(excitation);
            // throw std::invalid_argument("Voltage has not been processed");
        }
        else if (!excitation.get_voltage()->get_processed()->get_peak()) {
            process_voltage(excitation);
        }
        maximumVoltage = std::max(maximumVoltage, excitation.get_voltage()->get_processed()->get_peak().value());
    }

    return maximumVoltage;
}

double Inputs::get_maximum_voltage_rms(size_t windingIndex) {
    if (get_operating_points().size() == 0)
        throw std::invalid_argument("There are no operating points");

    double maximumVoltage = 0;
    for (auto& operatingPoint : get_operating_points()) {
        if (operatingPoint.get_excitations_per_winding().size() == 0)
            throw std::invalid_argument("There are no winding excitation in operating point");

        auto excitation = operatingPoint.get_excitations_per_winding()[windingIndex];
        if (!excitation.get_voltage()) 
            throw std::invalid_argument("Missing voltage in excitation");
        if (!excitation.get_voltage()->get_processed()) {
            process_voltage(excitation);
            // throw std::invalid_argument("Voltage has not been processed");
        }
        else if (!excitation.get_voltage()->get_processed()->get_rms()) {
            process_voltage(excitation);
        }
        maximumVoltage = std::max(maximumVoltage, excitation.get_voltage()->get_processed()->get_rms().value());
    }

    return maximumVoltage;
}

double Inputs::get_maximum_current_rms(size_t windingIndex) {
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

double Inputs::get_maximum_current_peak(size_t windingIndex) {
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

double Inputs::get_maximum_current_effective_frequency() {
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

double Inputs::get_maximum_current_dc_bias(size_t windingIndex) {
    if (get_operating_points().size() == 0)
        throw std::invalid_argument("There are no operating points");

    double maximumDcBias = 0;
    for (auto& operatingPoint : get_operating_points()) {
        if (operatingPoint.get_excitations_per_winding().size() == 0)
            throw std::invalid_argument("There are no winding excitation in operating point");

        auto excitation = operatingPoint.get_excitations_per_winding()[windingIndex];
        if (!excitation.get_current())
            throw std::invalid_argument("Missing current in excitation");
        if (!excitation.get_current()->get_processed())
            throw std::invalid_argument("Current has not been processed");
        maximumDcBias = std::max(maximumDcBias, excitation.get_current()->get_processed()->get_offset());
    }
    return maximumDcBias;
}

double Inputs::get_maximum_current_effective_frequency(size_t windingIndex) {
    if (get_operating_points().size() == 0)
        throw std::invalid_argument("There are no operating points");

    double maximumCurrentEffectiveFrequency = 0;
    for (auto& operatingPoint : get_operating_points()) {
        if (operatingPoint.get_excitations_per_winding().size() == 0)
            throw std::invalid_argument("There are no winding excitation in operating point");

        auto excitation = operatingPoint.get_excitations_per_winding()[windingIndex];
        if (!excitation.get_current()) 
            throw std::invalid_argument("Missing current in excitation");
        if (!excitation.get_current()->get_processed()) 
            throw std::invalid_argument("Current has not been processed");
        maximumCurrentEffectiveFrequency = std::max(maximumCurrentEffectiveFrequency, excitation.get_current()->get_processed()->get_effective_frequency().value());
    }
    return maximumCurrentEffectiveFrequency;
}

double Inputs::get_maximum_frequency() {
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

double Inputs::get_maximum_temperature() const {
    if (get_operating_points().size() == 0)
        throw std::invalid_argument("There are no operating points");

    double maximumTemperature = 0;
    for (const auto& operatingPoint : get_operating_points()) {
        maximumTemperature = std::max(maximumTemperature, operatingPoint.get_conditions().get_ambient_temperature());
    }

    return maximumTemperature;
}

DimensionWithTolerance Inputs::get_altitude() {
    if (!get_design_requirements().get_insulation())
        throw std::invalid_argument("Missing insulation in designRequirements");
    if (!get_design_requirements().get_insulation()->get_altitude())
        throw std::invalid_argument("Missing altitude in insulation requirements");

    return get_design_requirements().get_insulation()->get_altitude().value();
}

Cti Inputs::get_cti() {
    if (!get_design_requirements().get_insulation())
        throw std::invalid_argument("Missing insulation in designRequirements");
    if (!get_design_requirements().get_insulation()->get_cti())
        throw std::invalid_argument("Missing cti in insulation requirements");

    return get_design_requirements().get_insulation()->get_cti().value();
}

IsolationClass Inputs::get_insulation_type() {
    if (!get_design_requirements().get_insulation())
        throw std::invalid_argument("Missing insulation in designRequirements");
    if (!get_design_requirements().get_insulation()->get_insulation_type())
        throw std::invalid_argument("Missing insulation_type in insulation requirements");

    return get_design_requirements().get_insulation()->get_insulation_type().value();
}

DimensionWithTolerance Inputs::get_main_supply_voltage() {
    if (!get_design_requirements().get_insulation())
        throw std::invalid_argument("Missing insulation in designRequirements");
    if (!get_design_requirements().get_insulation()->get_main_supply_voltage())
        throw std::invalid_argument("Missing main_supply_voltage in insulation requirements");

    return get_design_requirements().get_insulation()->get_main_supply_voltage().value();
}

OvervoltageCategory Inputs::get_overvoltage_category() {
    if (!get_design_requirements().get_insulation())
        throw std::invalid_argument("Missing insulation in designRequirements");
    if (!get_design_requirements().get_insulation()->get_overvoltage_category())
        throw std::invalid_argument("Missing overvoltage_category in insulation requirements");

    return get_design_requirements().get_insulation()->get_overvoltage_category().value();
}

PollutionDegree Inputs::get_pollution_degree() {
    if (!get_design_requirements().get_insulation())
        throw std::invalid_argument("Missing insulation in designRequirements");
    if (!get_design_requirements().get_insulation()->get_pollution_degree())
        throw std::invalid_argument("Missing pollution_degree in insulation requirements");

    return get_design_requirements().get_insulation()->get_pollution_degree().value();
}

std::vector<InsulationStandards> Inputs::get_standards() {
    if (!get_design_requirements().get_insulation())
        throw std::invalid_argument("Missing insulation in designRequirements");
    if (!get_design_requirements().get_insulation()->get_standards()) {
        // ABT #1329: optional in MAS, but this operation cannot run without it. Never assume one.
        throw InvalidInputException(ErrorCode::MISSING_DATA,
            "designRequirements.insulation.standards is missing: insulation coordination (clearance, creepage, "
            "distance through insulation, lead sleeves) takes its tables from the standards the design names, and "
            "this design names none");
    }

    if (get_design_requirements().get_insulation()->get_standards()->empty()) {
        // ABT #1228: an empty list is no more a standard than an absent one; reading it as "no
        // standard asks for anything" returned 0 mm clearance and creepage.
        throw InvalidInputException(ErrorCode::MISSING_DATA,
            "designRequirements.insulation.standards is empty: insulation coordination takes its tables from the "
            "standards the design names, and this design names none");
    }
    return get_design_requirements().get_insulation()->get_standards().value();
}

bool Inputs::has_insulation_coordination_requirements() const {
    auto insulation = get_design_requirements().get_insulation();
    return insulation && insulation->get_standards() && !insulation->get_standards()->empty();
}

void Inputs::throw_if_json_lacks_required_fields(const json& inputsJson) {
    for (const char* field : {"designRequirements", "operatingPoints"}) {
        if (!inputsJson.is_object() || !inputsJson.contains(field)) {
            throw InvalidInputException(ErrorCode::MISSING_DATA,
                std::string("MAS inputs.") + field + " is missing; the MAS schema requires it");
        }
    }
}

WiringTechnology Inputs::get_wiring_technology() const {
    if (!get_design_requirements().get_wiring_technology()) {
        return defaults.wiringTechnology;
    }
    else {
        return get_design_requirements().get_wiring_technology().value();
    }
}

void Inputs::set_current_as_magnetizing_current(OperatingPoint* operatingPoint) {
    OperatingPointExcitation excitation = Inputs::get_primary_excitation(*operatingPoint);
    auto current = excitation.get_current().value();

    auto currentExcitation = excitation.get_current().value();
    auto currentExcitationWaveform = currentExcitation.get_waveform().value();
    if (!currentExcitation.get_processed() || !currentExcitation.get_harmonics()) {
        if (excitation.get_frequency() <= 0) {
            throw std::invalid_argument("Frequency has to be positive");
        }
        auto sampledCurrentWaveform = Inputs::calculate_sampled_waveform(currentExcitationWaveform, excitation.get_frequency());

        if (sampledCurrentWaveform.get_data().size() > 0 && ((sampledCurrentWaveform.get_data().size() & (sampledCurrentWaveform.get_data().size() - 1)) != 0)) {
            throw std::invalid_argument("sampledCurrentWaveform vector size is not a power of 2");
        }

        currentExcitation.set_harmonics(Inputs::calculate_harmonics_data(currentExcitationWaveform, sampledCurrentWaveform, excitation.get_frequency()));
        currentExcitation.set_processed(Inputs::calculate_processed_data(currentExcitation, sampledCurrentWaveform, true));
        excitation.set_current(currentExcitation);
    }
    excitation.set_magnetizing_current(excitation.get_current().value());
    operatingPoint->get_mutable_excitations_per_winding()[0] = excitation;
}

// Valleys of a record, found with a hysteresis of `noise`: a reversal smaller than that is not a
// flux reversal. Plateaus keep their first sample.
std::vector<size_t> Inputs::find_valleys(const std::vector<double>& data, double noise) {
    std::vector<size_t> valleys;
    int direction = 0;  // +1 rising towards a peak, -1 falling towards a valley, 0 not yet known
    size_t extremum = 0;
    for (size_t i = 1; i < data.size(); ++i) {
        if (direction == 0) {
            if (data[i] > data[0] + noise) {
                direction = 1;
                extremum = i;
            }
            else if (data[i] < data[0] - noise) {
                direction = -1;
                extremum = i;
            }
        }
        else if (direction > 0) {
            if (data[i] > data[extremum]) {
                extremum = i;
            }
            else if (data[extremum] - data[i] > noise) {
                direction = -1;
                extremum = i;
            }
        }
        else {
            if (data[i] < data[extremum]) {
                extremum = i;
            }
            else if (data[i] - data[extremum] > noise) {
                valleys.push_back(extremum);
                direction = 1;
                extremum = i;
            }
        }
    }
    return valleys;
}

std::optional<std::vector<size_t>> Inputs::find_switching_cycle_valleys(const std::vector<double>& time, const std::vector<double>& data) {
    if (data.size() != time.size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "Switching cycles: the record has " + std::to_string(data.size()) + " samples and " + std::to_string(time.size()) + " time points");
    }
    if (data.size() < 2 * minimumSwitchingCyclesInRecord + 1) {
        return std::nullopt;  // too few samples to hold that many cycles
    }
    auto [minimumIt, maximumIt] = std::minmax_element(data.begin(), data.end());
    double recordPeakToPeak = *maximumIt - *minimumIt;
    if (!(recordPeakToPeak > 0)) {
        return std::nullopt;
    }
    // A reversal below one part in 1e9 of the record's swing is the round-off of the waveform's
    // arithmetic (resampling, scaling current into flux), not a reversal.
    auto valleys = find_valleys(data, recordPeakToPeak * 1e-9);
    if (valleys.size() < minimumSwitchingCyclesInRecord + 1) {
        return std::nullopt;
    }
    double recordDuration = time.back() - time.front();
    for (size_t cycleIndex = 0; cycleIndex + 1 < valleys.size(); ++cycleIndex) {
        double duration = time[valleys[cycleIndex + 1]] - time[valleys[cycleIndex]];
        if (duration * static_cast<double>(minimumSwitchingCyclesInRecord) > recordDuration) {
            return std::nullopt;
        }
    }
    return valleys;
}

double Inputs::get_switching_frequency(const OperatingPointExcitation& excitation) {
    // The switching frequency of an excitation is its frequency, unless its current is a line-cycle
    // record: a record of repeated switching cycles (find_switching_cycle_valleys) described at the
    // frequency of its period, the line. Its switching frequency is then the rate of its cycles. The
    // record is read from its waveform at any excitation frequency, whether or not a processed copy
    // cached the waveform's spectrum (ABT #1764: a fixed f < 400 Hz gate gave a DMC record swept to 500 Hz
    // no line bias, and an unprocessed 50 Hz copy returned 50 Hz for want of cached harmonics: mu 300
    // instead of ~104 on an MPP 300 toroid, B 1.02 T instead of 0.35 T, 16.3 W instead of 1.74 W).
    double frequency = excitation.get_frequency();
    if (!excitation.get_current() || !excitation.get_current()->get_waveform()) {
        return frequency;
    }
    auto waveform = excitation.get_current()->get_waveform().value();
    if (!waveform.get_time()) {
        return frequency;  // an equidistant waveform is one period of its frequency, not a record
    }
    auto time = waveform.get_time().value();
    if (time.size() < 2) {
        return frequency;
    }
    // A record that spans that many periods of the excitation's frequency already carries the switching
    // frequency as its frequency (Kirchhoff's PFC/Vienna records do).
    if ((time.back() - time.front()) * frequency >= static_cast<double>(minimumSwitchingCyclesInRecord)) {
        return frequency;
    }
    auto valleys = find_switching_cycle_valleys(time, waveform.get_data());
    if (!valleys) {
        return frequency;
    }
    // The record's switching rate is its number of cycles over their span. A cycle in which the line
    // outran the ripple, so that it shows no valley, appears as one cycle twice the usual length: each
    // cycle counts as its length in typical cycles (the median one), rounded. A periodic record switches
    // at a harmonic of its period, so the rate is the nearest one: the cycles are bounded by samples, so
    // their span is off by up to a sample at each end.
    std::vector<double> durations;
    for (size_t cycleIndex = 0; cycleIndex + 1 < valleys->size(); ++cycleIndex) {
        durations.push_back(time[(*valleys)[cycleIndex + 1]] - time[(*valleys)[cycleIndex]]);
    }
    auto sortedDurations = durations;
    std::nth_element(sortedDurations.begin(), sortedDurations.begin() + sortedDurations.size() / 2, sortedDurations.end());
    double typicalDuration = sortedDurations[sortedDurations.size() / 2];
    double numberCycles = 0;
    for (auto duration : durations) {
        numberCycles += std::max(1.0, std::round(duration / typicalDuration));
    }
    double cyclesSpan = time[valleys->back()] - time[valleys->front()];
    double cycleRate = numberCycles / cyclesSpan;
    return std::round(cycleRate / frequency) * frequency;
}

double Inputs::get_magnetic_flux_density_peak(OperatingPointExcitation excitation, double switchingFrequency) {
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();

    if (!excitation.get_magnetic_flux_density()->get_harmonics() && magneticFluxDensity.get_waveform()) {
        auto magneticFluxDensityWaveform = magneticFluxDensity.get_waveform().value();
        auto sampledWaveform = Inputs::calculate_sampled_waveform(magneticFluxDensityWaveform, excitation.get_frequency());
        // Harmonic k at k times the excitation frequency, as every other MKF harmonics list (ABT #1880).
        magneticFluxDensity.set_harmonics(calculate_harmonics_data(magneticFluxDensityWaveform, sampledWaveform, excitation.get_frequency()));
        excitation.set_magnetic_flux_density(magneticFluxDensity);
    }

    if (excitation.get_magnetic_flux_density()->get_harmonics()) {
        auto harmonics = excitation.get_magnetic_flux_density()->get_harmonics().value();
        if (harmonics.get_amplitudes().size() > 1 && harmonics.get_frequencies()[1] != switchingFrequency) {
            for (size_t harmonicIndex = 1; harmonicIndex < harmonics.get_amplitudes().size(); ++harmonicIndex) {
                if (harmonics.get_frequencies()[harmonicIndex] == switchingFrequency) {
                    return harmonics.get_amplitudes()[harmonicIndex];
                }
            }
        }
    }

    return magneticFluxDensity.get_processed().value().get_peak().value();
}

double Inputs::get_magnetic_flux_density_peak_to_peak(OperatingPointExcitation excitation, double switchingFrequency) {
    auto magneticFluxDensity = excitation.get_magnetic_flux_density().value();

    if (!excitation.get_magnetic_flux_density()->get_harmonics() && magneticFluxDensity.get_waveform()) {
        auto magneticFluxDensityWaveform = magneticFluxDensity.get_waveform().value();
        auto sampledWaveform = Inputs::calculate_sampled_waveform(magneticFluxDensityWaveform, excitation.get_frequency());
        // Harmonic k at k times the excitation frequency, as every other MKF harmonics list (ABT #1880).
        magneticFluxDensity.set_harmonics(calculate_harmonics_data(magneticFluxDensityWaveform, sampledWaveform, excitation.get_frequency()));
        excitation.set_magnetic_flux_density(magneticFluxDensity);
    }

    if (excitation.get_magnetic_flux_density()->get_harmonics()) {
        auto harmonics = excitation.get_magnetic_flux_density()->get_harmonics().value();
        if (harmonics.get_amplitudes().size() > 1 && harmonics.get_frequencies()[1] != switchingFrequency) {
            for (size_t harmonicIndex = 1; harmonicIndex < harmonics.get_amplitudes().size(); ++harmonicIndex) {
                if (harmonics.get_frequencies()[harmonicIndex] == switchingFrequency) {
                    return harmonics.get_amplitudes()[harmonicIndex] * 2;
                }
            }
        }
    }

    return magneticFluxDensity.get_processed().value().get_peak_to_peak().value();
}

SignalDescriptor Inputs::get_current_with_effective_maximum(size_t windingIndex) {
    SignalDescriptor maximumCurrent;
    double maximumCurrentRmsTimesRootSquaredEffectiveFrequency = 0;
    for (size_t operatingPointIndex = 0; operatingPointIndex < get_operating_points().size(); ++operatingPointIndex) {
        if (!get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current()) {
            throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Current is missing");
        }
        if (!get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current()->get_processed()) {
            throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Current is not processed");
        }
        if (!get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current()->get_processed()->get_effective_frequency()) {
            auto current = get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current().value();
            if (current.get_waveform()) {
                auto processed = calculate_processed_data(current.get_waveform().value(), get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_frequency());
                current.set_processed(processed);
                get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex].set_current(current);
            }
            else {
                throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Current is not processed");
            }
        }
        double effectiveFrequency = get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current()->get_processed()->get_effective_frequency().value();
        double rms = get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current()->get_processed()->get_rms().value();

        auto currentRmsTimesRootSquaredEffectiveFrequency = rms * sqrt(effectiveFrequency);
        if (currentRmsTimesRootSquaredEffectiveFrequency > maximumCurrentRmsTimesRootSquaredEffectiveFrequency) {
            maximumCurrentRmsTimesRootSquaredEffectiveFrequency = currentRmsTimesRootSquaredEffectiveFrequency;
            maximumCurrent = get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current().value();
        }
    }

    return maximumCurrent;
}

std::vector<IsolationSide> Inputs::get_isolation_sides_used() {
    size_t numWindings = get_design_requirements().get_turns_ratios().size() + 1;
    
    if (!get_design_requirements().get_isolation_sides()) {
        // No isolation sides specified - auto-generate for all windings
        std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY};
        for (size_t windingIndex = 1; windingIndex < numWindings; ++windingIndex) {
            isolationSides.push_back(get_isolation_side_from_index(windingIndex));
        }
        get_mutable_design_requirements().set_isolation_sides(isolationSides);
    }
    else {
        // Isolation sides are provided - validate and extend if necessary
        auto isolationSides = get_design_requirements().get_isolation_sides().value();
        if (isolationSides.size() < numWindings) {
            // Extend isolation sides array for remaining windings
            // Assign them to SECONDARY side by default (most common case for multi-secondary transformers)
            for (size_t windingIndex = isolationSides.size(); windingIndex < numWindings; ++windingIndex) {
                isolationSides.push_back(get_isolation_side_from_index(windingIndex));
            }
            get_mutable_design_requirements().set_isolation_sides(isolationSides);
        }
    }

    std::vector<IsolationSide> isolationSidesUsed;
    std::vector<IsolationSide> isolationSidesFromRequirements = get_design_requirements().get_isolation_sides().value();
    for (auto isolationSide : isolationSidesFromRequirements) {
        if (std::find(isolationSidesUsed.begin(), isolationSidesUsed.end(), isolationSide) == isolationSidesUsed.end()) {
            isolationSidesUsed.push_back(isolationSide);
        }
    }

    return isolationSidesUsed;
}

} // namespace OpenMagnetics
