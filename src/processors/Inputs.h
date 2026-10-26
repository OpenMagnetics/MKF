#pragma once

#include "json.hpp"
#include "constructive_models/MasMigration.h"

#include <MAS.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <streambuf>
#include <vector>
#include "Definitions.h"
#include "support/Exceptions.h"
using json = nlohmann::json;

using namespace MAS;

namespace OpenMagnetics {

class Coil;

class Inputs : public MAS::Inputs {
  public:
    Inputs(json j, bool processWaveform = true, std::optional<std::variant<double, std::vector<double>>> magnetizingInductance = std::nullopt) {
        OpenMagnetics::compat::migrate_pre_1_0(j);
        throw_if_json_lacks_required_fields(j);
        from_json(j, *this);
        auto check_passed = check_integrity();
        if (!check_passed.first) {
            throw std::runtime_error("Missing inputs");
        }
        if (processWaveform) {
            process(magnetizingInductance);
        }
    }
    Inputs(std::vector<OperatingPoint> operatingPoints, DesignRequirements designRequirements, bool processWaveform = true, std::optional<std::variant<double, std::vector<double>>> magnetizingInductance = std::nullopt) {
        set_operating_points(operatingPoints);
        set_design_requirements(designRequirements);
        auto check_passed = check_integrity();
        if (!check_passed.first) {
            throw std::runtime_error("Missing inputs: " + check_passed.second);
        }
        if (processWaveform) {
            process(magnetizingInductance);
        }
    }
    Inputs() = default;
    virtual ~Inputs() = default;

    // ABT #1329: MAS requires inputs.designRequirements and inputs.operatingPoints. Without this
    // check a missing one surfaced as a raw "[json.exception.out_of_range.403] key
    // 'designRequirements' not found" from the generated parser, naming no MAS field.
    static void throw_if_json_lacks_required_fields(const json& inputsJson);
    // ABT #1329: insulation coordination (clearance, creepage, distance through insulation, lead
    // sleeves) is defined by the standards the design names. designRequirements.insulation.standards
    // is optional in MAS, so an insulation block without it is valid and simply names no standard
    // to coordinate by: true only when the insulation block carries a standards list. (An EMPTY
    // list, which MKF's own quick inputs build, keeps its previous meaning: coordinate by no table.)
    bool has_insulation_coordination_requirements() const;
    std::pair<bool, std::string> check_integrity();
    void process(std::optional<std::variant<double, std::vector<double>>> magnetizingInductance = std::nullopt);
    // isolationSides: the windings' sides; without them MKF's index convention applies (winding k on
    // get_isolation_side_from_index(k)). Windings on different sides oppose, so identical currents
    // on them cancel: that is a transformer, never a common-mode choke (ABT #1467).
    static OperatingPoint process_operating_point(OperatingPoint operatingPoint, double magnetizingInductance, std::optional<std::vector<double>> turnsRatios = std::nullopt, bool isDmcTopology = false,
                                                  std::optional<std::vector<IsolationSide>> isolationSides = std::nullopt);

    static bool is_waveform_sampled(const Waveform& waveform);
    static bool is_waveform_imported(const Waveform& waveform);
    static bool is_multiport_inductor(const OperatingPoint& operatingPoint, std::optional<std::vector<IsolationSide>> isolationSides = std::nullopt);
    // isolationSides, when given: a common-mode choke has every winding on one side. Windings on
    // different sides oppose, so identical currents on them cancel (a transformer), ABT #1467.
    static bool can_be_common_mode_choke(const OperatingPoint& operatingPoint, std::optional<std::vector<IsolationSide>> isolationSides = std::nullopt);
    static bool can_be_differential_mode_choke(const OperatingPoint& operatingPoint, bool isDmcTopology = false);

    static Waveform calculate_sampled_waveform(Waveform waveform, double frequency, std::optional<size_t> numberPoints=std::nullopt, std::optional<size_t> maximumNumberPoints=std::nullopt);
    static ProcessedWaveform calculate_processed_data(Waveform waveform, double frequency, bool includeAdvancedData=true, std::optional<ProcessedWaveform> processed=std::nullopt);
    static ProcessedWaveform calculate_processed_data(SignalDescriptor excitation, Waveform sampledWaveform, bool includeAdvancedData=true, std::optional<ProcessedWaveform> processed=std::nullopt);
    static ProcessedWaveform calculate_processed_data(Harmonics harmonics, Waveform waveform, bool includeAdvancedData=true, std::optional<ProcessedWaveform> processed=std::nullopt);
    static Harmonics calculate_harmonics_data(const Waveform& waveform, double frequency);
    // Harmonics of `waveform` from its own knots when it has a time axis (the exact series, no
    // sampling), else of `sampledWaveform`, its uniform resample: a data-only waveform has no
    // knots to be exact about (ABT #1460).
    static Harmonics calculate_harmonics_data(const Waveform& waveform, const Waveform& sampledWaveform, double frequency);
    static OperatingPointExcitation prune_harmonics(OperatingPointExcitation excitation, double windingLossesHarmonicAmplitudeThreshold, std::optional<size_t> mainHarmonicIndex=std::nullopt);
    static SignalDescriptor prune_harmonics(SignalDescriptor signalDescriptor, double windingLossesHarmonicAmplitudeThreshold, std::optional<size_t> mainHarmonicIndex=std::nullopt);
    static OperatingPoint prune_harmonics(OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold, std::optional<size_t> mainHarmonicIndex=std::nullopt);

    // static OperatingPointExcitation reflect_waveforms(OperatingPointExcitation excitation, double ratio);
    static SignalDescriptor reflect_waveform(SignalDescriptor excitation, double ratio);
    static SignalDescriptor reflect_waveform(SignalDescriptor signal, double ratio, WaveformLabel label);
    // Secondary excitation of an ideal two-winding transformer with turns ratio n = N1/N2:
    // voltage v1/n (Faraday, any shape), current reflected by its label. The result carries
    // secondaryName (or no name), never the primary's. Throws if the primary lacks either signal.
    static OperatingPointExcitation calculate_reflected_secondary(OperatingPointExcitation primaryExcitation, double turnRatio, std::optional<std::string> secondaryName = std::nullopt);
    // Primary excitation of the same transformer, from the secondary's: voltage v2*n (Faraday,
    // any shape), current reflected by its label with ratio 1/n, so a flyback secondary current
    // comes back as a flyback primary one. The result carries primaryName (or no name), never
    // the secondary's. Throws if the secondary lacks either signal.
    static OperatingPointExcitation calculate_reflected_primary(OperatingPointExcitation secondaryExcitation, double turnRatio, std::optional<std::string> primaryName = std::nullopt);

    static bool is_standardized(const SignalDescriptor& signal);
    static SignalDescriptor standardize_waveform(SignalDescriptor parameter, double frequency);
    static Waveform reconstruct_signal(Harmonics harmonics, double frequency);
    OperatingPoint get_operating_point(size_t index);
    const OperatingPoint& get_operating_point_ref(size_t index) const { return get_operating_points()[index]; }
    OperatingPointExcitation get_winding_excitation(size_t operatingPointIndex=0, size_t windingIndex=0);
    OperatingPointExcitation get_primary_excitation(size_t operatingPointIndex=0);
    // Perf (Phase 6): const-ref overload avoids deep-copying the OperatingPoint
    // just to read a single excitation. Hot in AreaProduct and Loss filters.
    static OperatingPointExcitation get_primary_excitation(const OperatingPoint& operatingPoint);

    static SignalDescriptor calculate_induced_voltage(OperatingPointExcitation& excitation, double magnetizingInductance, bool compress=true);

    /**
     * @brief Peak flux linkage of one period of the voltage waveform.
     *
     * Returns (max − min)/2 of the running ∫(V − ⟨V⟩) dt over one period:
     * the peak of the flux linkage about its mean, NOT its peak-to-peak
     * (a sine of amplitude V_pk gives V_pk/ω, as the fallback). This is the V·s
     * quantity that pairs with Faraday's law (`B_peak = V·s / (N · A_e)`)
     * for non-sinusoidal voltages — DAB squares, CLLLC trapezoids,
     * push-pull rectangulars, etc. — where the sinusoidal V_peak/ω
     * shortcut underestimates B_peak by ~36 %.
     *
     * Fallback: if the operating point doesn't carry a sampled voltage
     * waveform but has a processed peak + frequency, return V_peak/ω
     * (the sinusoidal estimate). Returns 0 if neither is available.
     *
     * Iterates over ALL operating points and ALL excitations (or the
     * single OP / primary excitation overload below), returning the
     * worst-case excursion.
     */
    static double calculate_max_volt_seconds(const OperatingPointExcitation& excitation);
    static double calculate_max_volt_seconds(const OperatingPoint& operatingPoint);
    // Every winding's current sampled on ONE time grid over the period (numberPoints, or the
    // finest power of 2 the currents need), and whether it carries current at each instant
    // (|i| above 1 % of that winding's peak). std::nullopt when no winding has a current (a voltage-only
    // operating point); throws when some windings have one and another has none (ABT #1814).
    struct WindingConduction {
        std::vector<std::vector<double>> currents;
        std::vector<std::vector<bool>> conducting;
    };
    static std::optional<WindingConduction> calculate_winding_conduction(const OperatingPoint& operatingPoint, std::optional<size_t> numberPoints = std::nullopt);
    // The instants at which a winding's current is the magnetizing current: the primary conducting
    // while no other winding does, any other winding conducting while the primary does not.
    static std::vector<bool> calculate_magnetizing_conduction_instants(const WindingConduction& conduction, size_t windingIndex);
    // Whether the windings' currents give the magnetizing DC: a single winding, or windings that do
    // not always conduct together (flyback, forward reset). Time-aware (ABT #1467).
    static bool include_dc_offset_into_magnetizing_current(const OperatingPoint& operatingPoint, const std::vector<double>& turnsRatios);
    // The DC of a winding's magnetizing current, taken where that winding carries the magnetizing
    // MMF alone (ABT #907, #1467); 0 when it never does. Throws when a winding has no current (ABT #1814).
    static double calculate_magnetizing_current_dc_anchor(const OperatingPoint& operatingPoint, size_t windingIndex, Waveform voltageSampledWaveform, double magnetizingInductance);
    static double calculate_magnetizing_current_dc_anchor(OperatingPointExcitation excitation, const WindingConduction& conduction, size_t windingIndex,
                                                          Waveform voltageSampledWaveform, double magnetizingInductance);
    // The ampere-turn current of an operating point: i_m(t) = sum_k c_k (N_k/N_ref) i_k(t), the
    // single current that, in the reference winding, drives the same core flux as all windings
    // together. Built from the winding current WAVEFORMS on one common period grid, so phase and
    // dot direction count: c_k = +1 for a primary-side winding, -1 otherwise (MAS's excitation
    // sign convention: primary-side currents passive, the others source), N_k/N_ref from the
    // design's turnsRatios (Np/Nk against winding 0, re-referred to the reference) or 1:1 when
    // none are given; the reference is the first primary-side winding. Its processed values come
    // from the summed waveform (offset = signed DC, peak = magnitude peak). Throws when a winding
    // has no current waveform (processed values alone carry no phase), when the windings'
    // frequencies differ, or when the isolation sides are not given for every winding.
    static SignalDescriptor calculate_ampere_turn_current(const OperatingPoint& operatingPoint,
                                                          const std::vector<double>& turnsRatios,
                                                          const std::vector<IsolationSide>& isolationSides);

    // Largest net volt-second imbalance over one period, as a fraction of the AC flux-linkage
    // swing, that a voltage may carry and still be integrated into a magnetizing current (ABT #1585).
    // A periodic steady-state winding voltage has zero net volt-seconds; what remains in a simulated
    // or measured waveform is sampling error (at most 0.12 of the swing on the multi-winding and
    // switching-frequency single-winding cases of the Henry corpus, 2026-10-02). The integral of an
    // unbalanced voltage drifts by the imbalance every period, so above this bound the integrated
    // current describes the drift, not the magnetizing current: the line-cycle PFC inductor
    // voltages that motivated it carried 16-46 times their swing and gave kA peaks.
    static constexpr double maximumVoltSecondImbalance = 0.25;
    // The volt-second balance of a sampled voltage over its period (ABT #1855): its mean, the net
    // volt-seconds that mean leaves per period, the AC flux-linkage swing of the rest, and the net
    // as a fraction of that swing. A constant voltage has no swing (acSwing 0, fraction not finite
    // when its mean is not zero).
    struct VoltSecondBalance {
        double meanVoltage;
        double period;
        double netVoltSeconds;
        double acSwing;
        double imbalance;
    };
    static VoltSecondBalance calculate_volt_second_balance(const Waveform& voltageSampledWaveform, const std::string& context);
    // The per-topic imbalance check (ABT #1855), shared by check_volt_second_balance and by
    // ResultsAnalyzer: MKF's message when the voltage's net volt-seconds over its period exceed
    // maximumVoltSecondImbalance of its AC flux-linkage swing (or it is a non-zero constant), nullopt
    // when it is balanced. `context` names the operating point and winding for the message.
    static std::optional<std::string> describe_volt_second_imbalance(const Waveform& voltageSampledWaveform, const std::string& context);
    // Throws InvalidInputException with describe_volt_second_imbalance's message when the voltage is
    // not balanced. Since ABT #1855 only the calculation of an induced current asked for that voltage
    // (calculate_magnetizing_current of an excitation, the bindings' calculate_induced_current) and
    // the advisers' up-front input check refuse it; loading, autocompleting and simulating a design
    // integrate the voltage's balanced part and ResultsAnalyzer reports the imbalance.
    static void check_volt_second_balance(const Waveform& voltageSampledWaveform, const std::string& context);
    // ABT #2002: a winding's terminal voltage is its own resistive drop plus what moves the flux,
    // v(t) = R_dc i(t) + N dPhi/dt. A DC current through the winding's resistance gives the voltage a
    // mean that is not flux drift, so the volt-second balance finding judges the voltage with that
    // drop taken out. A magnetizing current is derived from the voltage as given, which MKF treats as
    // the induced voltage (converter models and the web give ideal voltages). `remove_resistive_drop`
    // returns the excitation with its voltage replaced by v(t) - R_dc i(t), exact on the knots of
    // both waveforms, over one period (the voltage's own sampled axis, or [t0, t0 + T]); it throws when
    // the excitation has no voltage or no current, or `dcResistance` is not a finite value >= 0.
    static OperatingPointExcitation remove_resistive_drop(const OperatingPointExcitation& excitation, double dcResistance);
    // Each winding's DC resistance at `temperature` from its wire and turns (WindingOhmicLosses::
    // calculate_dc_resistance_per_winding), or nullopt when the coil has not been wound: without turns
    // MKF cannot know the winding's length, hence its resistance. A wound coil whose resistance cannot
    // be computed throws.
    static std::optional<std::vector<double>> get_winding_dc_resistances(const Coil& coil, double temperature);
    // The mean of the voltage AS GIVEN (not its resampled copy). A winding in periodic steady state
    // has zero net volt-seconds, so what remains (probe offset, I*R drop, a capture that is not a
    // whole number of periods, or a voltage that is not balanced at all) is subtracted before
    // integrating the magnetizing current, and logged as a warning above
    // voltSecondImbalanceWarningThreshold. It does not refuse an unbalanced voltage (ABT #1855).
    static constexpr double voltSecondImbalanceWarningThreshold = 0.01;
    static double calculate_steady_state_mean_voltage(const OperatingPointExcitation& excitation,
                                                      const Waveform& voltageSampledWaveform,
                                                      const std::string& context);
    // ABT #1718: the advisers' loss ranking (MagnetizingInductance::calculate_inductance_and_magnetic_flux_density)
    // integrates the primary voltage of every operating point into a magnetizing current, and that
    // integration rejects an unbalanced voltage (check_volt_second_balance). Reached only after every
    // candidate has been wound, the rejection killed the whole adviser late, or culled every candidate.
    // This runs the same check once, up front, on exactly the excitations that path integrates, and
    // throws an InvalidInputException starting "Adviser inputs are invalid:" that names the operating
    // point and the winding. Excitations routed to a path that does not integrate the voltage (a
    // single-winding current, a common-mode choke, a multiport inductor, a flyback/unipolar
    // triangular current, an open-core family model) are not checked. What the adviser does not know
    // yet is assumed in the direction that checks less: `coilIsolationSides` nullopt = no coil yet,
    // `coreShapeFamily` / `coreType` nullopt = no core yet.
    static void check_adviser_inputs_volt_second_balance(const Inputs& inputs,
                                                         std::optional<std::vector<IsolationSide>> coilIsolationSides = std::nullopt,
                                                         std::optional<std::string> primaryWindingName = std::nullopt,
                                                         std::optional<CoreShapeFamily> coreShapeFamily = std::nullopt,
                                                         std::optional<CoreType> coreType = std::nullopt);

    // A single-winding inductor's winding current IS its magnetizing current (ABT #1585); the
    // voltage is integrated only when the current is absent or several windings share the flux.
    // Returns the current with complete processed data (peak included), or nullopt when the rule
    // does not apply. Every place that fills magnetizingCurrent must go through this first.
    static std::optional<SignalDescriptor> get_single_winding_magnetizing_current(const OperatingPointExcitation& excitation,
                                                                                  size_t numberWindings);

    // A winding's magnetizing current, the single decision every caller shares (ABT #1467):
    // Faraday on its voltage (voltageSampledWaveform, or its own voltage sampled), DC from the
    // currents where include_dc_offset_into_magnetizing_current says they give it.
    static SignalDescriptor calculate_magnetizing_current(const OperatingPoint& operatingPoint,
                                                            size_t windingIndex,
                                                            double magnetizingInductance,
                                                            const std::vector<double>& turnsRatios,
                                                            bool compress,
                                                            std::optional<Waveform> voltageSampledWaveform = std::nullopt);
    static SignalDescriptor calculate_magnetizing_current(OperatingPointExcitation& excitation,
                                                            double magnetizingInductance,
                                                            bool compress,
                                                            double dcCurrent=0);
    static SignalDescriptor calculate_magnetizing_current(OperatingPointExcitation& excitation,
                                                            Waveform voltageSampledWaveform,
                                                            double magnetizingInductance,
                                                            bool compress,
                                                            double dcCurrent=0);
    // The two excitation overloads above are the induced current asked for one voltage (the
    // bindings' calculate_induced_current): they refuse a voltage that is not volt-second balanced
    // (check_volt_second_balance), whose message the web shows next to it. Everything that loads,
    // autocompletes or simulates a design integrates through these instead, which take the voltage's
    // balanced part (its mean removed) and do not refuse it (ABT #1855).
    static SignalDescriptor integrate_magnetizing_current(OperatingPointExcitation& excitation,
                                                          double magnetizingInductance,
                                                          bool compress,
                                                          double dcCurrent);
    static SignalDescriptor integrate_magnetizing_current(OperatingPointExcitation& excitation,
                                                          Waveform voltageSampledWaveform,
                                                          double magnetizingInductance,
                                                          bool compress,
                                                          double dcCurrent);
    // A bool in the dcCurrent slot is the removed (addOffset, alternatingConduction) signature: it
    // converted silently to a 1 A DC (PyOpenMagnetics did exactly that). Callers that want the
    // currents' DC use the OperatingPoint overload above.
    static SignalDescriptor calculate_magnetizing_current(OperatingPointExcitation& excitation, double magnetizingInductance, bool compress, bool dcCurrent) = delete;
    static SignalDescriptor calculate_magnetizing_current(OperatingPointExcitation& excitation, Waveform voltageSampledWaveform, double magnetizingInductance, bool compress, bool dcCurrent) = delete;
    static SignalDescriptor integrate_magnetizing_current(OperatingPointExcitation& excitation, double magnetizingInductance, bool compress, bool dcCurrent) = delete;
    static SignalDescriptor integrate_magnetizing_current(OperatingPointExcitation& excitation, Waveform voltageSampledWaveform, double magnetizingInductance, bool compress, bool dcCurrent) = delete;
    static SignalDescriptor add_offset_to_excitation(SignalDescriptor signalDescriptor,
                                                             double offset,
                                                             double frequency);
    static OperatingPointExcitation get_excitation_with_proportional_current(OperatingPointExcitation excitation, double proportion);
    static OperatingPointExcitation get_excitation_with_proportional_voltage(OperatingPointExcitation excitation, double proportion);

    static void make_waveform_size_power_of_two(OperatingPoint* operatingPoint);

    static Inputs create_quick_operating_point(double frequency,
                                                      double magnetizingInductance,
                                                      double temperature,
                                                      WaveformLabel waveShape,
                                                      double peakToPeak,
                                                      double dutyCycle,
                                                      double dcCurrent,
                                                      std::vector<double> turnsRatios={});
    static Inputs create_quick_operating_point_only_current(double frequency,
                                                                  double magnetizingInductance,
                                                                  double temperature,
                                                                  WaveformLabel waveShape,
                                                                  double peakToPeak,
                                                                  double dutyCycle,
                                                                  double dcCurrent,
                                                                  std::vector<double> turnsRatios={});
    static Inputs create_quick_operating_point_only_current(double frequency,
                                                                  double magnetizingInductance,
                                                                  double temperature,
                                                                  WaveformLabel waveShape,
                                                                  std::vector<double> peakToPeaks,
                                                                  double dutyCycle,
                                                                  double dcCurrent,
                                                                  std::vector<double> turnsRatios={});
    static OperatingPoint create_operating_point_with_sinusoidal_current_mask(double frequency,
                                                                              double magnetizingInductance,
                                                                              double temperature,
                                                                              std::vector<double> turnsRatios,
                                                                              std::vector<double> currentPeakMask,
                                                                              double currentOffset = 0);

    static WaveformLabel try_guess_waveform_label(const Waveform& waveform);
    // ABT #1330: whether a processed description defines a waveform on its own. A signal given only
    // by its processed parameters (label, peakToPeak, offset, dutyCycle...) is rebuilt with
    // create_waveform, which reads the duty cycle for every non-sinusoidal label and the dead time
    // for the *_WITH_DEADTIME ones; CUSTOM and RECTANGULAR_DCM have no parametric form at all.
    // Throws naming what is missing (`what` says which signal), never assumes a value.
    static void throw_if_processed_cannot_define_waveform(const ProcessedWaveform& processed, const std::string& what);
    static Waveform create_waveform(ProcessedWaveform processed, double frequency);
    // phase: radians, applied to SINUSOIDAL waveforms only (positive phase
    // advances the waveform — e.g. phase=pi/2 turns sin(wt) into sin(wt+pi/2)=cos(wt),
    // so the resulting signal leads the same-frequency sine by 90 degrees).
    // Ignored for non-sinusoidal labels (use `skew` for those).
    static Waveform create_waveform(WaveformLabel label, double peakToPeak, double frequency, double dutyCycle=0.5, double offset=0, double deadTime=0, double skew=0, double phase=0);
    // givenWaveform: the waveform as its producer gave it, when `waveform` is a sampled or compressed copy (ABT #1878).
    static ProcessedWaveform calculate_basic_processed_data(const Waveform& waveform, const std::optional<Waveform>& givenWaveform = std::nullopt);
    // See WaveformProcessor::compress_waveform: with the frequency a closed axis is not closed twice.
    static Waveform compress_waveform(const Waveform& waveform, std::optional<double> frequency = std::nullopt);

    static Waveform calculate_derivative_waveform(Waveform waveform);
    static Waveform calculate_integral_waveform(Waveform waveform, bool subtractAverage=false);
    // Time average over the waveform's period. Without a frequency the axis is taken as closed
    // ([t0, t0 + T]) and integrated by trapezoids. With the excitation frequency an axis that is
    // open by one step ([t0, t0 + T - dt], the sampled convention, uniform only) is recognised
    // (WaveformProcessor::check_time_axis_spans_one_period) and averaged as the sampled waveform it
    // is: each sample holds until the next, the last until t0 + T, so the mean of N uniform samples
    // is their sum / N, the same convention calculate_integral_waveform integrates with (ABT #1713).
    static double calculate_waveform_average(const Waveform& waveform, std::optional<double> frequency = std::nullopt);
    static Waveform multiply_waveform(Waveform waveform, double scalarValue);
    static Waveform sum_waveform(Waveform waveform, double scalarValue);
    // Pointwise add two waveforms. Both must have the same number of data points.
    // Time vector is taken from `a` (must be aligned with `b`).
    static Waveform add_waveforms(Waveform a, const Waveform& b);
    // Pointwise subtract: returns a - b. Same length requirement.
    static Waveform subtract_waveforms(Waveform a, const Waveform& b);
    // Circularly shift waveform samples by `samples` positions to the right
    // (i.e. delay by samples * dt). Negative values shift left (advance).
    // Operates on the data vector only; time vector (if any) is preserved.
    static Waveform shift_waveform(Waveform waveform, int samples);
    // SignalDescriptor wrappers: perform the op on the underlying sampled
    // waveform and recompute harmonics and processed data at `frequency`.
    static SignalDescriptor add_signals(const SignalDescriptor& a, const SignalDescriptor& b, double frequency);
    static SignalDescriptor subtract_signals(const SignalDescriptor& a, const SignalDescriptor& b, double frequency);
    static SignalDescriptor shift_signal(const SignalDescriptor& s, int samples, double frequency);

    static double try_guess_duty_cycle(const Waveform& waveform, WaveformLabel label=WaveformLabel::CUSTOM, double frequency=0);
    static double calculate_instantaneous_power(OperatingPointExcitation excitation);

    static double calculate_waveform_coefficient(OperatingPoint* operatingPoint);

    static void scale_time_to_frequency(Inputs& inputs, double newFrequency, bool cleanFrequencyDependentFields=false, bool processSignals=false, bool useCurrentAsBase=true);
    static void scale_time_to_frequency(OperatingPoint& operatingPoint, double newFrequency, bool cleanFrequencyDependentFields=false, bool processSignals=false, bool useCurrentAsBase=true);
    static void scale_time_to_frequency(OperatingPointExcitation& excitation, double newFrequency, bool cleanFrequencyDependentFields=false, bool processSignals=false, bool useCurrentAsBase=true);
    static Waveform scale_time_to_frequency(Waveform waveform, double newFrequency);
    
    void set_operating_point_by_index(const OperatingPoint& value, size_t index) {
        get_mutable_operating_points()[index] = value;
    }
    static void set_current_as_magnetizing_current(OperatingPoint* operatingPoint);

    // A record needs at least this many complete reversal cycles, none of them longer than this
    // many-th part of the record, to be taken as repeated switching cycles (ABT #1636): a line-cycle
    // record. A periodic excitation is one period of its frequency: one loop, or a few reversals
    // (ringing, multi-level or burst waveforms, Gibbs ripple of a reconstructed step), at least one of
    // which spans a large part of the period.
    static constexpr size_t minimumSwitchingCyclesInRecord = 16;
    // Valleys of a record, found with a hysteresis of `noise`. Plateaus keep their first sample.
    static std::vector<size_t> find_valleys(const std::vector<double>& data, double noise);
    // The valleys that bound a record's switching cycles, or nothing when the record is not a record
    // of repeated switching cycles by the definition above.
    static std::optional<std::vector<size_t>> find_switching_cycle_valleys(const std::vector<double>& time, const std::vector<double>& data);
    static double get_switching_frequency(const OperatingPointExcitation& excitation);
    static double get_magnetic_flux_density_peak(OperatingPointExcitation excitation, double switchingFrequency);
    static double get_magnetic_flux_density_peak_to_peak(OperatingPointExcitation excitation, double switchingFrequency);
    // A multiport inductor's magnetizing current (ABT #1467). Flyback-labelled (the windings conduct
    // alternately): Faraday on the primary voltage, DC anchored at the primary current's peak
    // (ABT #907). Without a primary voltage, or with every winding on one isolation side: the
    // ampere-turn sum of the winding currents, referred to the primary (exact while they conduct).
    static SignalDescriptor get_multiport_inductor_magnetizing_current(const OperatingPoint& operatingPoint,
                                                                       double magnetizingInductance,
                                                                       const std::vector<double>& turnsRatios,
                                                                       std::optional<std::vector<IsolationSide>> isolationSides);
    // Whether get_multiport_inductor_magnetizing_current takes the magnetizing current from the primary voltage
    // (Faraday, a flyback with a primary voltage) rather than from the windings' ampere-turn sum (ABT #1942).
    static bool is_multiport_inductor_magnetizing_current_from_voltage(const OperatingPointExcitation& primaryExcitation);
    // The primary current carries a FLYBACK_PRIMARY / FLYBACK_SECONDARY label (processed or ancillary).
    static bool has_flyback_label(const OperatingPointExcitation& excitation);
    static SignalDescriptor get_common_mode_choke_magnetizing_current(const OperatingPoint& operatingPoint);
    static SignalDescriptor get_differential_mode_choke_magnetizing_current(const OperatingPoint& operatingPoint);

    DimensionWithTolerance get_altitude();
    Cti get_cti();
    IsolationClass get_insulation_type();
    DimensionWithTolerance get_main_supply_voltage();
    OvervoltageCategory get_overvoltage_category();
    PollutionDegree get_pollution_degree();
    std::vector<InsulationStandards> get_standards();
    WiringTechnology get_wiring_technology() const;
    double get_maximum_voltage_peak();
    double get_maximum_voltage_rms();
    double get_maximum_current_effective_frequency();
    double get_maximum_current_peak();
    double get_maximum_current_rms();
    double get_maximum_frequency();
    double get_maximum_temperature() const;
    double get_maximum_voltage_peak(size_t windingIndex);
    double get_maximum_voltage_rms(size_t windingIndex);
    double get_maximum_current_peak(size_t windingIndex);
    double get_maximum_current_rms(size_t windingIndex);
    double get_maximum_current_effective_frequency(size_t windingIndex);
    double get_maximum_current_dc_bias(size_t windingIndex = 0);
    SignalDescriptor get_current_with_effective_maximum(size_t windingIndex = 0);
    std::vector<IsolationSide> get_isolation_sides_used();

};
void from_json(const json& j, Inputs& x);
void to_json(json& j, const Inputs& x);
void to_file(std::filesystem::path filepath, const Inputs & x);

inline void from_json(const json& j, Inputs& x) {
    x.set_design_requirements(j.at("designRequirements").get<DesignRequirements>());
    x.set_operating_points(j.at("operatingPoints").get<std::vector<OperatingPoint>>());
}
inline void to_json(json& j, const Inputs& x) {
    j = json::object();
    j["designRequirements"] = x.get_design_requirements();
    j["operatingPoints"] = x.get_operating_points();
}
inline void to_file(std::filesystem::path filepath, const Inputs & x) {
    json masJson;
    to_json(masJson, x);

    std::ofstream myfile;
    myfile.open(filepath);
    myfile << masJson;
}


// What a magnetic does in its circuit, for the choices that are categorical rather than a magnitude:
// whether the core is gapped, whether the inductor loss optimiser applies, whether an ungappable
// ferrite toroid is dropped, and the common-mode / differential-mode choke handling. The
// magnitudes (saturation, the turn seed, the energy check, the loss-sweep cap) never read it:
// they use the magnetizing current and the volt-seconds of every design.
enum class MagneticRole { STORING, TRANSFORMER, COMMON_MODE_CHOKE, DIFFERENTIAL_MODE_CHOKE };

// Which evidence decided the role, most reliable first. SINGLE_WINDING: one winding stores energy.
// TOPOLOGY: designRequirements.topology, when the topology table knows it (a common- or
// differential-mode choke topology names the component, so it decides even before SINGLE_WINDING). APPLICATION:
// designRequirements.subApplication (common-/differential-mode noise filtering, power filtering).
// ISOLATION_SIDES: every declared winding on one isolation side. STORAGE_RATIO: the energy the
// windings exchange with the core within a period against the energy that passes through it,
// from the voltage and current of every winding. NET_MAGNETOMOTIVE_FORCE: the net ampere-turns
// against the largest winding's, from the currents and turns ratios alone. NO_CURRENTS: voltages
// only, or not every winding's current: a transformer.
enum class MagneticRoleTier { SINGLE_WINDING, TOPOLOGY, APPLICATION, ISOLATION_SIDES, STORAGE_RATIO, NET_MAGNETOMOTIVE_FORCE, NO_CURRENTS };

struct MagneticRoleDecision {
    MagneticRole role;
    MagneticRoleTier tier;
    // The storage ratio (kappa), the largest over the operating points, when the STORAGE_RATIO tier ran.
    std::optional<double> storageRatio;
    // The net magnetomotive force ratio (rho), the largest over the operating points, when the
    // NET_MAGNETOMOTIVE_FORCE tier ran.
    std::optional<double> netMagnetomotiveForceRatio;
    // What the decision rests on and anything the caller should know (e.g. the winding signs read
    // in the passive convention because the declared sides did not balance the power).
    std::vector<std::string> findings;
};

// A storage ratio or net magnetomotive force ratio at or above this is an energy-storing magnetic.
inline constexpr double magneticRoleStoringThreshold = 0.5;

// One winding and no turns ratios.
bool is_single_winding_design(const Inputs& inputs);
// The role a converter topology gives its magnetic; nullopt for a topology the table does not know.
std::optional<MagneticRole> topology_magnetic_role(MAS::Topology topology);
// The role from the declarations alone (tiers SINGLE_WINDING to ISOLATION_SIDES); nullopt when
// none decides. The only tiers that can yield a common-mode or differential-mode choke.
std::optional<MagneticRoleDecision> declared_magnetic_role(const Inputs& inputs);
// The role, through every tier. Never throws for missing descriptive metadata.
MagneticRoleDecision magnetic_role(const Inputs& inputs);
// The topology the magnetizing-inductance model is given: the declared one, else a
// differential-mode choke when the declarations make the magnetic one (its magnetizing current is
// the sum of its windings' line currents).
std::optional<MAS::Topology> magnetizing_model_topology(const Inputs& inputs);

bool operator==(Inputs lhs, Inputs rhs);

inline bool operator==(Inputs lhs, Inputs rhs) {
    // TODO: Add more comparisond
    bool isEqual = resolve_dimensional_values(lhs.get_design_requirements().get_magnetizing_inductance()) == resolve_dimensional_values(rhs.get_design_requirements().get_magnetizing_inductance()) &&
                   lhs.get_operating_points().size() == rhs.get_operating_points().size();

    if (isEqual) {
        // FIXED: BUG-09 — Check sizes match
        if (lhs.get_design_requirements().get_turns_ratios().size() != rhs.get_design_requirements().get_turns_ratios().size()) return false;
        for (size_t i = 0; i < lhs.get_design_requirements().get_turns_ratios().size(); ++i) {
            isEqual &= resolve_dimensional_values(lhs.get_design_requirements().get_turns_ratios()[i]) == resolve_dimensional_values(rhs.get_design_requirements().get_turns_ratios()[i]);
        }
        for (size_t i = 0; i < lhs.get_operating_points().size(); ++i) {
            isEqual &= lhs.get_operating_points()[i].get_excitations_per_winding().size() == rhs.get_operating_points()[i].get_excitations_per_winding().size();
        }
    }

    return isEqual;
}

} // namespace OpenMagnetics
