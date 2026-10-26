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

    // Throws std::invalid_argument unless the axis spans exactly one period of `frequency` (within
    // 0.1 %), closed [t0, t0 + T] or, on a uniform axis, open [t0, t0 + T - dt]; true when open.
    static bool check_time_axis_spans_one_period(const std::vector<double>& time, double frequency);

    static MAS::Waveform calculate_sampled_waveform(MAS::Waveform waveform,
                                                    double frequency,
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

    // Complex Fourier coefficients c_0..c_K of one period of a waveform, by harmonic index k
    // (frequency k * excitation frequency), referred to the waveform's first instant:
    //   x(t) = c_0 + sum_k 2 Re(c_k exp(+j 2 pi k (t - timeOrigin) / period)).
    // The ONE source of every harmonic's amplitude (calculate_harmonics_data: |c_0|, 2 |c_k|) and
    // phase (calculate_harmonic_phasor), so the two can never disagree (ABT #1880). The method
    // and period rule are those listed under calculate_harmonics_data below: open uniform axis
    // = samples (period N dt), closed axis = knots (period = span), denser = resampled, no axis
    // = plain DFT of the samples (period 1 / frequency).
    struct FourierSeries {
        std::vector<std::complex<double>> coefficients;
        double timeOrigin;  // s, the instant the coefficients are referred to
        double period;      // s, the period the coefficients were computed over
    };
    static FourierSeries calculate_fourier_series(const MAS::Waveform& waveform,
                                                  double frequency,
                                                  size_t numberPointsSampledWaveforms = 128);

    // Peak phasor X of the harmonic at harmonicFrequency, x(t) = Re(X exp(+j w t)) with t the
    // waveform's own (absolute) time: X = 2 c_k (X = c_0 at DC) from calculate_fourier_series, k =
    // harmonicFrequency / frequency, rotated from the series' time origin to t = 0. |X| is the
    // MAS amplitude of that harmonic. Throws when harmonicFrequency is not a multiple of frequency
    // (to 1e-6 of the fundamental) or k is beyond the series. A waveform without a time axis is
    // sampled over one period first, as Inputs::calculate_harmonics_data(waveform, sampled, f)
    // does when it lists its harmonics.
    static std::complex<double> calculate_harmonic_phasor(const FourierSeries& series,
                                                          double frequency,
                                                          double harmonicFrequency);
    static std::complex<double> calculate_harmonic_phasor(const MAS::Waveform& waveform,
                                                          double frequency,
                                                          double harmonicFrequency,
                                                          size_t numberPointsSampledWaveforms = 128);
    // The series calculate_harmonic_phasor(waveform, ...) takes its phasors from: a waveform
    // without a time axis is sampled over one period first.
    static FourierSeries calculate_waveform_fourier_series(const MAS::Waveform& waveform,
                                                           double frequency,
                                                           size_t numberPointsSampledWaveforms = 128);

    // Relative tolerance a MAS harmonics list and its waveform are held to (ABT #1880): a listed
    // amplitude a needs a waveform phasor |X| >= harmonicConsistencyTolerance * a to have a phase.
    static constexpr double harmonicConsistencyTolerance = 1e-6;
    // Largest MAS amplitude of the series, max(|c_0|, 2 |c_k|): the scale a harmonic is zero on.
    static double calculate_largest_harmonic_amplitude(const FourierSeries& series);
    // A harmonic that is zero both in the list (listedAmplitude) and in the waveform (phasor):
    // both below harmonicConsistencyTolerance times the waveform's largest harmonic. It carries
    // no current, so it contributes nothing whatever its phase, and it has no phase to take.
    static bool is_zero_harmonic(const FourierSeries& series, double listedAmplitude, std::complex<double> phasor);

    // Harmonics 0..K of a waveform (ABT #1460):
    //  - with a time axis and fewer knots than numberPointsSampledWaveforms: the exact series of
    //    its piecewise-linear knots (calculate_exact_fourier_coefficients), K = N/2;
    //  - a sampled waveform (is_waveform_sampled: N uniform samples on [t0, t0 + N dt)): the
    //    exact series of the periodic linear interpolant of the samples, which is the N-point
    //    DFT times sinc^2(pi k / N), K = N/2; N points on an uneven time axis, or on a CLOSED
    //    uniform axis ([t0, t0 + T], last point = first), are knots, as above (ABT #1879);
    //  - denser knots than N that are not a sampled waveform: sampled at their own power-of-2
    //    resolution first (a closed form over thousands of knots costs O(knots x harmonics)),
    //    then as above;
    //  - no time axis: samples only, no knots to be exact about: the plain DFT, harmonics
    //    0..N/2-1, as before.
    // The coefficients come from calculate_fourier_series, the same ones that give the phases
    // (calculate_harmonic_phasor): amplitude and phase of a harmonic are one complex number.
    static MAS::Harmonics calculate_harmonics_data(const MAS::Waveform& waveform,
                                                   double frequency,
                                                   bool trimHarmonics = true,
                                                   double harmonicAmplitudeThreshold = 0.05,
                                                   size_t numberPointsSampledWaveforms = 128);

    static MAS::ProcessedWaveform calculate_processed_data(MAS::Waveform waveform,
                                                           double frequency,
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
                                                           size_t numberPointsSampledWaveforms = 128,
                                                           std::optional<MAS::Waveform> givenWaveform = std::nullopt);
    // `waveform` is what is measured (peaks, the label guess); `givenWaveform` is the waveform as
    // its producer gave it, before any sampling, when that differs. A label's vertex layout is
    // checked on, and its offset and duty are read off, the given waveform when it is not a
    // sample; a sampled or compressed copy is never held to the label's layout (ABT #1878).
    static MAS::ProcessedWaveform calculate_basic_processed_data(const MAS::Waveform& waveform,
                                                                 size_t numberPointsSampledWaveforms = 128,
                                                                 const std::optional<MAS::Waveform>& givenWaveform = std::nullopt);

    // Waveform-shape helpers (ported alongside the DSP entry points above).
    static bool is_waveform_sampled(const MAS::Waveform& waveform, size_t numberPointsSampledWaveforms = 128);
    static bool is_waveform_imported(const MAS::Waveform& waveform, size_t numberPointsSampledWaveforms = 128);
    static double calculate_waveform_average(const MAS::Waveform& waveform);
    // Keeps only the corners of a waveform (points where the slope changes) and closes it
    // periodically. Without a frequency the axis is taken as sampled, open by one step
    // ([t0, t0 + T - dt]), and the first sample is appended at t_N-1 + dt. With the excitation
    // frequency a closed axis ([t0, t0 + T], which may also hold 128 points) is recognised: it
    // already ends at t0 + T and gets no extra point; appending one stretched it to T + dt
    // (ABT #1467). Any other axis is closed as before; the span is not validated here.
    static MAS::Waveform compress_waveform(const MAS::Waveform& waveform, std::optional<double> frequency = std::nullopt);
    static MAS::WaveformLabel try_guess_waveform_label(const MAS::Waveform& waveform, size_t numberPointsSampledWaveforms = 128);
    static double try_guess_duty_cycle(const MAS::Waveform& waveform,
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
