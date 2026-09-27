#pragma once

#include <MAS.hpp>

#include <complex>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace OpenMagnetics {

// WaveformProcessor: the pure converter-waveform DSP extracted out of
// OpenMagnetics::Inputs / OpenMagnetics::Topology. It depends ONLY on the MAS
// data types (<MAS.hpp>) and the standard library, with NO coupling to the
// magnetics stack (Magnetic, MagnetizingInductance, Painter) and NO dependency
// on the MKF Settings singleton. Every value that used to be read from
// Settings is now an explicit function parameter, defaulted to the same value
// the singleton is constructed with (numberPointsSampledWaveforms = 128,
// trimHarmonics = true, harmonicAmplitudeThreshold = 0.05). This lets the class
// be reused by the Kirchhoff converter-model library without dragging in MKF.
//
// Behavior is identical to the original Inputs/Topology implementations; the
// Inputs and Topology methods now delegate here, passing settings.get_...().
class WaveformProcessor {
  public:
    static MAS::Waveform create_waveform(MAS::ProcessedWaveform processed,
                                         double frequency,
                                         size_t numberPointsSampledWaveforms = 128);

    // phase: radians, applied to SINUSOIDAL waveforms only (positive phase
    // advances the waveform). Ignored for non-sinusoidal labels (use `skew`).
    static MAS::Waveform create_waveform(MAS::WaveformLabel label,
                                         double peakToPeak,
                                         double frequency,
                                         double dutyCycle = 0.5,
                                         double offset = 0,
                                         double deadTime = 0,
                                         double skew = 0,
                                         double phase = 0,
                                         size_t numberPointsSampledWaveforms = 128);

    static MAS::Waveform calculate_sampled_waveform(MAS::Waveform waveform,
                                                    double frequency = 0,
                                                    std::optional<size_t> numberPoints = std::nullopt,
                                                    size_t numberPointsSampledWaveforms = 128,
                                                    std::optional<size_t> maximumNumberPoints = std::nullopt);

    // Exact complex Fourier coefficients c_0..c_numberHarmonics of the PERIODIC piecewise-linear
    // waveform through the knots (time[i], data[i]), period = time.back() - time.front(), the
    // period closing with a step from data.back() to data.front() when they differ.
    // Closed form per segment, no sampling, so no aliasing (ABT #1460):
    //   c_k = 1/(j 2 pi k) * sum_segments dv * exp(-j w_k t_mid) * sinc(w_k dt / 2)
    // A repeated sample (dt = 0, dv = 0) contributes nothing; a step (dt = 0, dv != 0) is the
    // dt -> 0 limit of the same term, exactly. Time running backwards throws.
    // Amplitudes in MAS convention are |c_0| and 2 |c_k|.
    static std::vector<std::complex<double>> calculate_exact_fourier_coefficients(const std::vector<double>& time,
                                                                                  const std::vector<double>& data,
                                                                                  size_t numberHarmonics);

    // Harmonics 0..K of a waveform (ABT #1460):
    //  - with a time axis and fewer knots than numberPointsSampledWaveforms: the exact series of
    //    its piecewise-linear knots (calculate_exact_fourier_coefficients), K = N/2;
    //  - a sampled waveform (is_waveform_sampled: N uniform samples on [t0, t0 + N dt)): the
    //    exact series of the periodic linear interpolant of the samples, which is the N-point
    //    DFT times sinc^2(pi k / N), K = N/2; N points on an uneven time axis are knots, as above;
    //  - denser knots than N that are not a sampled waveform: sampled at their own power-of-2
    //    resolution first (a closed form over thousands of knots costs O(knots x harmonics)),
    //    then as above;
    //  - no time axis: samples only, no knots to be exact about: the plain DFT, harmonics
    //    0..N/2-1, as before.
    static MAS::Harmonics calculate_harmonics_data(MAS::Waveform waveform,
                                                   double frequency,
                                                   bool trimHarmonics = true,
                                                   double harmonicAmplitudeThreshold = 0.05,
                                                   size_t numberPointsSampledWaveforms = 128);

    static MAS::ProcessedWaveform calculate_processed_data(MAS::Waveform waveform,
                                                           std::optional<double> frequency = std::nullopt,
                                                           bool includeAdvancedData = true,
                                                           std::optional<MAS::ProcessedWaveform> processed = std::nullopt,
                                                           bool trimHarmonics = true,
                                                           double harmonicAmplitudeThreshold = 0.05,
                                                           size_t numberPointsSampledWaveforms = 128);
    static MAS::ProcessedWaveform calculate_processed_data(MAS::SignalDescriptor excitation,
                                                           MAS::Waveform sampledWaveform,
                                                           bool includeAdvancedData = true,
                                                           std::optional<MAS::ProcessedWaveform> processed = std::nullopt,
                                                           size_t numberPointsSampledWaveforms = 128);
    static MAS::ProcessedWaveform calculate_processed_data(MAS::Harmonics harmonics,
                                                           MAS::Waveform waveform,
                                                           bool includeAdvancedData = true,
                                                           std::optional<MAS::ProcessedWaveform> processed = std::nullopt,
                                                           size_t numberPointsSampledWaveforms = 128);
    static MAS::ProcessedWaveform calculate_basic_processed_data(MAS::Waveform waveform,
                                                                 size_t numberPointsSampledWaveforms = 128);

    // Waveform-shape helpers (ported alongside the DSP entry points above).
    static bool is_waveform_sampled(MAS::Waveform waveform, size_t numberPointsSampledWaveforms = 128);
    static bool is_waveform_imported(MAS::Waveform waveform, size_t numberPointsSampledWaveforms = 128);
    static double calculate_waveform_average(MAS::Waveform waveform);
    static MAS::Waveform compress_waveform(const MAS::Waveform& waveform);
    static MAS::WaveformLabel try_guess_waveform_label(MAS::Waveform waveform, size_t numberPointsSampledWaveforms = 128);
    static double try_guess_duty_cycle(MAS::Waveform waveform,
                                       MAS::WaveformLabel label = MAS::WaveformLabel::CUSTOM,
                                       double frequency = 0,
                                       size_t numberPointsSampledWaveforms = 128);

    static MAS::OperatingPointExcitation complete_excitation(MAS::Waveform currentWaveform,
                                                             MAS::Waveform voltageWaveform,
                                                             double switchingFrequency,
                                                             std::string name = "",
                                                             bool trimHarmonics = true,
                                                             double harmonicAmplitudeThreshold = 0.05,
                                                             size_t numberPointsSampledWaveforms = 128);
};

} // namespace OpenMagnetics
