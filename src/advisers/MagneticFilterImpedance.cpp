#include "advisers/MagneticFilter.h"
#include "advisers/MagneticFilterInternal.h"
#include "constructive_models/NumberTurns.h"
#include "physical_models/ComplexPermeability.h"
#include "physical_models/Impedance.h"
#include "physical_models/WindingLosses.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "support/Exceptions.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <numbers>

namespace OpenMagnetics {

// Phase 5: extracted from MagneticFilter.cpp.
// This translation unit owns the impedance-family filters:
//   - MagneticFilterCoreMinimumImpedance (the heavyweight CMC suppression filter)
//   - MagneticFilterEffectiveResistance
//   - MagneticFilterProximityFactor
//   - MagneticFilterImpedance (post-cap evaluator)
// Shared helpers (prepare_bobbin_for_non_pqi etc.) live in
// advisers/MagneticFilterInternal.h.

std::pair<bool, double> MagneticFilterCoreMinimumImpedance::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    const auto& core = magnetic->get_core();

    double primaryCurrentRms = 0;
    // Phase 6 (perf): cache operating-points by const-ref to avoid OperatingPoint deep copies.
    const auto& operatingPoints = inputs->get_operating_points();
    for (size_t operatingPointIndex = 0; operatingPointIndex < operatingPoints.size(); ++operatingPointIndex) {
        primaryCurrentRms = std::max(primaryCurrentRms, Inputs::get_primary_excitation(operatingPoints[operatingPointIndex]).get_current().value().get_processed().value().get_rms().value());
    }

    std::string shapeName = core.get_shape_name();
    prepare_bobbin_for_non_pqi(magnetic, shapeName);

    // For impedance filter, start searching from 1 turn rather than from inductance-based
    // initial value. In interference suppression applications, impedance is the primary
    // requirement and we need to find the optimal turn count that meets impedance while
    // keeping SRF above the operating frequency range.
    NumberTurns numberTurns(1);

    Coil coil = magnetic->get_coil();

    double conductingArea = primaryCurrentRms / defaults.maximumCurrentDensity;
    auto wire = Wire::get_wire_for_conducting_area(conductingArea, defaults.ambientTemperature, false);
    coil.get_mutable_functional_description()[0].set_wire(wire);
    coil.unwind();

    if (!inputs->get_design_requirements().get_minimum_impedance()) {
        throw InvalidInputException("Minimum impedance missing from requirements");
    }

    auto minimumImpedanceRequirement = inputs->get_design_requirements().get_minimum_impedance().value();
    auto windingWindowArea = magnetic->get_mutable_coil().resolve_bobbin().get_winding_window_area();

    // ABT #121.2: the complex-permeability impedance solve is temperature
    // dependent, so evaluate it at the operating point's actual ambient (like
    // MagneticFilterTemperatureRise) instead of the hardcoded 25 degC default.
    // When there is genuinely no operating point (a pure-impedance CMC spec),
    // there is no ambient to read — keep the documented default rather than
    // throwing on a path that never carried a temperature.
    double ambientTemperature = inputs->get_operating_points().size() > 0
                                    ? inputs->get_maximum_temperature()
                                    : defaults.ambientTemperature;

    // Compute the peak magnetizing current for DC bias correction.
    // For multi-winding chokes (DMC/CMC), the magnetizing current is set
    // by the converter model as the vector sum of all winding currents.
    double magnetizingCurrentPeak = 0;
    double effectiveLength = core.get_processed_description()->get_effective_parameters().get_effective_length();
    if (inputs->get_operating_points().size() > 0) {
        auto& primaryExcitation = inputs->get_operating_points()[0].get_excitations_per_winding()[0];
        if (primaryExcitation.get_magnetizing_current() &&
            primaryExcitation.get_magnetizing_current()->get_processed() &&
            primaryExcitation.get_magnetizing_current()->get_processed()->get_peak()) {
            magnetizingCurrentPeak = primaryExcitation.get_magnetizing_current()->get_processed()->get_peak().value();
        } else if (primaryExcitation.get_current() &&
                   primaryExcitation.get_current()->get_processed() &&
                   primaryExcitation.get_current()->get_processed()->get_peak()) {
            magnetizingCurrentPeak = primaryExcitation.get_current()->get_processed()->get_peak().value();
        }
    }

    bool validDesign = true;
    bool validMaterial = true;
    double totalImpedanceExtra = 0;
    int timeout = defaults.coilAdviserCmcImpedanceMaxIterations;
    bool jumpDefinitive = false;
    int64_t jumpedRequiredN = 0;

    // Analytical first-jump: inductive impedance scales as N² (with mild
    // N-dependence via DC-bias and stray capacitance). For each impedance
    // requirement, compute Z(N=1) once, then jump straight to
    // N_required = ceil(sqrt(Z_target / Z(N=1))). For first-pass filtering
    // this is definitive: the do-loop below is skipped (perf) and we accept
    // the analytical answer. Capacitance/SRF rolloff is a second-order
    // effect that downstream filters re-validate on the survivor cap.
    {
        Coil unitCoil = coil;
        unitCoil.get_mutable_functional_description()[0].set_number_turns(1.0);
        int64_t requiredN = 1;
        double extraSum = 0;
        bool jumpOk = true;
        for (auto impedanceAtFrequency : minimumImpedanceRequirement) {
            auto frequency = impedanceAtFrequency.get_frequency();
            auto required = impedanceAtFrequency.get_impedance().get_magnitude();
            try {
                double zUnit;
                if (magnetizingCurrentPeak > 0 && effectiveLength > 0) {
                    double H_dc_unit = 1.0 * magnetizingCurrentPeak / effectiveLength;
                    zUnit = abs(_impedanceModel.calculate_impedance(core, unitCoil, frequency, H_dc_unit, ambientTemperature));
                } else {
                    zUnit = abs(_impedanceModel.calculate_impedance(core, unitCoil, frequency));
                }
                if (zUnit <= 0 || !std::isfinite(zUnit)) { jumpOk = false; break; }
                int64_t nNeeded = static_cast<int64_t>(std::ceil(std::sqrt(required / zUnit)));
                if (nNeeded < 1) nNeeded = 1;
                if (nNeeded > requiredN) requiredN = nNeeded;
                // Use the post-jump impedance estimate for the score
                double estZ = zUnit * static_cast<double>(requiredN) * static_cast<double>(requiredN);
                extraSum += std::max(0.0, estZ - required);
            } catch (const ModelNotAvailableException&) {
                jumpOk = false;
                break;
            }
        }
        if (jumpOk) {
            // Geometric feasibility check: do the turns physically fit?
            double wireOuterRadius = resolve_dimensional_values(wire.get_outer_diameter().value()) / 2.0;
            double conductorAreaTotal = static_cast<double>(requiredN) * std::numbers::pi * wireOuterRadius * wireOuterRadius;
            if (conductorAreaTotal >= windingWindowArea) {
                magnetic->set_coil(std::move(coil));
                return {false, 0};
            }
            coil.get_mutable_functional_description()[0].set_number_turns(static_cast<double>(requiredN));
            // Validation pass at the analytical N. The Z ∝ N² assumption used
            // for the jump above breaks for materials with field-dependent
            // permeability (powder cores: µ_r drops as H_dc = N·I_pk/lₑ
            // grows with N), so the analytical N often under-shoots Z_target
            // by 10–30 %. When that happens, do a Newton-style proportional
            // re-bump — N_new = ceil(N · sqrt(Z_target / Z_measured)) — and
            // re-validate, using the same _impedanceModel. Capped at a few
            // iterations because each iteration tightens the gap by the
            // d log Z / d log N slope (≈ 2 in the linear regime, ≈ 1 in deep
            // saturation), so 3–5 attempts converge for any physical core.
            constexpr int kMaxRebumpIterations = 5;
            double measuredExtra = 0;
            bool meetsAll = false;
            for (int iter = 0; iter < kMaxRebumpIterations; ++iter) {
                measuredExtra = 0;
                meetsAll = true;
                double worstShortfallRatio = 1.0;  // measured / required at the worst frequency
                bool modelMissing = false;
                for (auto impedanceAtFrequency : minimumImpedanceRequirement) {
                    auto frequency = impedanceAtFrequency.get_frequency();
                    auto required = impedanceAtFrequency.get_impedance().get_magnitude();
                    try {
                        double impedance;
                        if (magnetizingCurrentPeak > 0 && effectiveLength > 0) {
                            double H_dc = static_cast<double>(requiredN) * magnetizingCurrentPeak / effectiveLength;
                            impedance = abs(_impedanceModel.calculate_impedance(core, coil, frequency, H_dc, ambientTemperature));
                        } else {
                            impedance = abs(_impedanceModel.calculate_impedance(core, coil, frequency));
                        }
                        if (impedance < required) {
                            meetsAll = false;
                            double ratio = impedance / required;  // < 1 here
                            if (ratio < worstShortfallRatio) worstShortfallRatio = ratio;
                        } else {
                            measuredExtra += (impedance - required);
                        }
                    } catch (const ModelNotAvailableException&) {
                        modelMissing = true;
                        break;
                    }
                }
                if (modelMissing) {
                    meetsAll = false;
                    break;
                }
                if (meetsAll) break;
                // Proportional re-bump using the worst-frequency shortfall.
                // Multiplying N by sqrt(1/ratio) restores Z ∝ N² scaling
                // exactly when µ is constant; with field-dependent µ it
                // converges geometrically at the local d log Z / d log N
                // slope.
                int64_t bumpedN = static_cast<int64_t>(std::ceil(static_cast<double>(requiredN) * std::sqrt(1.0 / worstShortfallRatio)));
                if (bumpedN <= requiredN) bumpedN = requiredN + 1;  // ensure forward progress
                requiredN = bumpedN;
                // Re-check geometric feasibility before the next probe.
                double conductorAreaTotalIter = static_cast<double>(requiredN) * std::numbers::pi * wireOuterRadius * wireOuterRadius;
                if (conductorAreaTotalIter >= windingWindowArea) {
                    meetsAll = false;
                    break;
                }
                coil.get_mutable_functional_description()[0].set_number_turns(static_cast<double>(requiredN));
            }
            if (!meetsAll) {
                magnetic->set_coil(std::move(coil));
                return {false, 0};
            }
            jumpDefinitive = true;
            jumpedRequiredN = requiredN;
            totalImpedanceExtra = measuredExtra;
        }
    }

    if (!jumpDefinitive) do {
        totalImpedanceExtra = 0;
        validDesign = true;
        auto numberTurnsCombination = numberTurns.get_next_number_turns_combination();

        if (numberTurnsCombination[0] * std::numbers::pi * pow(resolve_dimensional_values(wire.get_outer_diameter().value()) / 2, 2) >= windingWindowArea) {

            validMaterial = false;
            break;
        }
        coil.get_mutable_functional_description()[0].set_number_turns(numberTurnsCombination[0]);
        auto selfResonantFrequency = _impedanceModel.calculate_self_resonant_frequency(core, coil);

        for (auto impedanceAtFrequency : minimumImpedanceRequirement) {
            auto frequency = impedanceAtFrequency.get_frequency();
            if (frequency > defaults.selfResonantFrequencyMargin * selfResonantFrequency) {

                validDesign = false;
                break;
            }
        }

        if (!validDesign) {
            break;
        }

        for (auto impedanceAtFrequency : minimumImpedanceRequirement) {
            auto frequency = impedanceAtFrequency.get_frequency();
            auto minimumImpedanceRequired = impedanceAtFrequency.get_impedance();
            try {
                double impedance;
                if (magnetizingCurrentPeak > 0 && effectiveLength > 0) {
                    double H_dc = numberTurnsCombination[0] * magnetizingCurrentPeak / effectiveLength;
                    impedance = abs(_impedanceModel.calculate_impedance(core, coil, frequency, H_dc, defaults.ambientTemperature));
                } else {
                    impedance = abs(_impedanceModel.calculate_impedance(core, coil, frequency));
                }

                if (impedance < minimumImpedanceRequired.get_magnitude()) {
                    validDesign = false;
                    break;
                }
                else {
                    totalImpedanceExtra += (impedance - minimumImpedanceRequired.get_magnitude());
                }

            }
            catch (const ModelNotAvailableException &exc) {

                validMaterial = false;
            }
        }

        timeout--;
    }
    while(!validDesign && validMaterial && timeout > 0);

    // Reference jumpedRequiredN to keep parity with the original code and
    // silence the unused-variable warning that lived here pre-extraction.
    (void)jumpedRequiredN;

    if (validDesign && validMaterial) {
        // Skip the expensive fast_wind() / are_sections_and_layers_fitting()
        // validation here. With ~thousands of candidate cores in suppression
        // flows this winding work dominates the runtime even though most
        // cores will be culled by the downstream cap. The cheap
        // winding-area-vs-conductor-area check inside the do-loop already
        // rejects geometrically infeasible designs; the final
        // correct_windings() pass after the cap performs full winding
        // validation on the surviving top-N.
        magnetic->set_coil(std::move(coil));
        return {true, totalImpedanceExtra};
    }

    magnetic->set_coil(std::move(coil));
    return {false, 0};
}

std::pair<bool, double> MagneticFilterEffectiveResistance::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;

    if (inputs->get_operating_points().size() > 0 && magnetic->get_mutable_coil().get_functional_description().size() != inputs->get_operating_points()[0].get_excitations_per_winding().size()) {
        return {false, 0.0};        
    }

    for (size_t windingIndex = 0; windingIndex < magnetic->get_coil().get_functional_description().size(); ++windingIndex) {
        auto winding = magnetic->get_coil().get_functional_description()[windingIndex];
        auto maximumEffectiveFrequency = inputs->get_maximum_current_effective_frequency(windingIndex);
        auto temperature = inputs->get_maximum_temperature();
        auto [auxValid, auxScoring] = evaluate_magnetic(winding, maximumEffectiveFrequency, temperature);
        valid &= auxValid;
        scoring += auxScoring;
    }
    scoring /= magnetic->get_coil().get_functional_description().size();

    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterEffectiveResistance::evaluate_magnetic(Winding winding, double effectivefrequency, double temperature) {
    auto wire = Coil::resolve_wire(winding);

    double effectiveResistancePerMeter = WindingLosses::calculate_effective_resistance_per_meter(wire, effectivefrequency, temperature);

    double valid = true;
    // double valid = effectiveResistancePerMeter < defaults.maximumEffectiveCurrentDensity;

    return {valid, effectiveResistancePerMeter};
}

std::pair<bool, double> MagneticFilterProximityFactor::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;

    if (inputs->get_operating_points().size() > 0 && magnetic->get_mutable_coil().get_functional_description().size() != inputs->get_operating_points()[0].get_excitations_per_winding().size()) {
        return {false, 0.0};        
    }

    for (size_t windingIndex = 0; windingIndex < magnetic->get_coil().get_functional_description().size(); ++windingIndex) {
        auto winding = magnetic->get_coil().get_functional_description()[windingIndex];
        auto maximumEffectiveFrequency = inputs->get_maximum_current_effective_frequency(windingIndex);
        auto temperature = inputs->get_maximum_temperature();
        double effectiveSkinDepth = WindingSkinEffectLosses::calculate_skin_depth("copper", maximumEffectiveFrequency, temperature);
        auto [auxValid, auxScoring] = evaluate_magnetic(winding, effectiveSkinDepth, temperature);
        valid &= auxValid;
        scoring += auxScoring;
    }
    scoring /= magnetic->get_coil().get_functional_description().size();

    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterProximityFactor::evaluate_magnetic(Winding winding, double effectiveSkinDepth, double temperature) {
    auto wire = Coil::resolve_wire(winding);

    if (!wire.get_number_conductors()) {
        wire.set_number_conductors(1);
    }
    double proximityFactor = wire.get_minimum_conducting_dimension() / effectiveSkinDepth * pow(wire.get_number_conductors().value() * winding.get_number_parallels() * winding.get_number_turns(), 2);
    // proximityFactor = wire.get_minimum_conducting_dimension() / effectiveSkinDepth * pow(winding.get_number_parallels() / (std::max(wire.get_maximum_outer_width(), wire.get_maximum_outer_height())), 2);

    double valid = true;
    // double valid = effectiveResistancePerMeter < defaults.maximumEffectiveCurrentDensity;

    return {valid, proximityFactor};
}

std::optional<MagneticFilterImpedance::MeasuredImpedanceCurve> MagneticFilterImpedance::get_measured_impedance_curve(const Magnetic& magnetic) {
    auto manufacturerInfo = magnetic.get_manufacturer_info();
    if (!manufacturerInfo || !manufacturerInfo->get_datasheet_info() || !manufacturerInfo->get_datasheet_info()->get_electrical()) {
        return std::nullopt;
    }
    std::string reference = manufacturerInfo->get_reference() ? manufacturerInfo->get_reference().value() : magnetic.get_reference();
    std::optional<MeasuredImpedanceCurve> curve;
    // Named copies throughout: these getters return optionals by value, and ranging over a
    // temporary's value() would iterate a destroyed object.
    const auto electricalEntries = manufacturerInfo->get_datasheet_info()->get_electrical().value();
    for (const auto& electrical : electricalEntries) {
        if (electrical.get_subtype() != ElectricalSubtype::COMMON_MODE_CHOKE || !electrical.get_impedance_points()) {
            continue;
        }
        MeasuredImpedanceCurve entryCurve;
        std::optional<std::optional<std::string>> curveWinding;
        const auto impedancePoints = electrical.get_impedance_points().value();
        for (const auto& point : impedancePoints) {
            if (point.get_current()) {
                double current = point.get_current().value();
                if (!std::isfinite(current)) {
                    throw InvalidDatasheetImpedanceException(reference, "a point at " + std::to_string(point.get_frequency()) +
                                                             " Hz has a non-finite DC-bias current");
                }
                if (current != 0) {
                    continue;  // a DC-biased measurement: not the zero-bias curve
                }
            }
            double frequency = point.get_frequency();
            double magnitude = point.get_impedance().get_magnitude();
            if (!std::isfinite(frequency) || frequency <= 0) {
                throw InvalidDatasheetImpedanceException(reference, "a point has frequency " + std::to_string(frequency) +
                                                         " Hz; a measured frequency must be finite and positive");
            }
            if (!std::isfinite(magnitude) || magnitude <= 0) {
                throw InvalidDatasheetImpedanceException(reference, "the point at " + std::to_string(frequency) + " Hz has |Z| " +
                                                         std::to_string(magnitude) + " Ohm; a measured magnitude must be finite and positive");
            }
            if (!curveWinding) {
                curveWinding = point.get_winding();
            }
            else if (curveWinding.value() != point.get_winding()) {
                throw InvalidDatasheetImpedanceException(reference, "the zero-bias points belong to more than one winding, so there is more than one curve to choose from");
            }
            entryCurve.push_back({frequency, magnitude});
        }
        if (entryCurve.empty()) {
            continue;
        }
        if (curve) {
            throw InvalidDatasheetImpedanceException(reference, "more than one commonModeChoke entry carries a zero-bias impedance curve");
        }
        std::sort(entryCurve.begin(), entryCurve.end());
        for (size_t index = 1; index < entryCurve.size(); ++index) {
            if (entryCurve[index].first == entryCurve[index - 1].first) {
                throw InvalidDatasheetImpedanceException(reference, "two zero-bias points are measured at " + std::to_string(entryCurve[index].first) + " Hz");
            }
        }
        curve = std::move(entryCurve);
    }
    return curve;
}

double MagneticFilterImpedance::interpolate_measured_impedance(const MeasuredImpedanceCurve& curve, double frequency) {
    if (curve.empty() || !(frequency >= curve.front().first && frequency <= curve.back().first)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Measured impedance asked at " + std::to_string(frequency) +
                                    " Hz, outside the measured range" +
                                    (curve.empty() ? std::string(" (no points)") : " " + std::to_string(curve.front().first) + "-" +
                                                                                    std::to_string(curve.back().first) + " Hz") +
                                    "; the measurement is never extrapolated");
    }
    auto upper = std::lower_bound(curve.begin(), curve.end(), frequency,
                                  [](const std::pair<double, double>& point, double value) { return point.first < value; });
    if (upper->first == frequency) {
        return upper->second;
    }
    auto lower = std::prev(upper);
    double fraction = (std::log(frequency) - std::log(lower->first)) / (std::log(upper->first) - std::log(lower->first));
    return std::exp(std::log(lower->second) + fraction * (std::log(upper->second) - std::log(lower->second)));
}

MagneticFilterImpedance::RequirementCoverage MagneticFilterImpedance::classify_requirement(Magnetic* magnetic, Inputs* inputs) {
    RequirementCoverage coverage;
    auto [minimumMaterialFrequency, maximumMaterialFrequency] = ComplexPermeability().get_frequency_range(magnetic->get_core().resolve_material());
    coverage.minimumMaterialFrequency = minimumMaterialFrequency;
    coverage.maximumMaterialFrequency = maximumMaterialFrequency;
    // A named copy: get_minimum_impedance() returns the optional by value, so ranging over its
    // value() directly would iterate a destroyed temporary.
    const auto requirement = inputs->get_design_requirements().get_minimum_impedance().value();
    bool measuredCurveRead = false;
    for (const auto& impedanceAtFrequency : requirement) {
        double frequency = impedanceAtFrequency.get_frequency();
        if (frequency >= minimumMaterialFrequency && frequency <= maximumMaterialFrequency) {
            coverage.sources.push_back(ImpedanceSource::MODEL);
            continue;
        }
        // Outside the material's data: the part's own measurement, read only when a point needs it.
        if (!measuredCurveRead) {
            coverage.measuredCurve = get_measured_impedance_curve(*magnetic);
            measuredCurveRead = true;
        }
        if (coverage.measuredCurve && frequency >= coverage.measuredCurve->front().first && frequency <= coverage.measuredCurve->back().first) {
            coverage.sources.push_back(ImpedanceSource::MEASURED);
        }
        else {
            coverage.sources.push_back(ImpedanceSource::NOT_JUDGED);
        }
    }
    return coverage;
}

std::optional<std::vector<double>> MagneticFilterImpedance::get_judged_frequencies(Magnetic* magnetic, Inputs* inputs) const {
    if (!inputs->get_design_requirements().get_minimum_impedance()) {
        return std::nullopt;
    }
    auto coverage = classify_requirement(magnetic, inputs);
    const auto requirement = inputs->get_design_requirements().get_minimum_impedance().value();
    std::vector<double> judgedFrequencies;
    for (size_t index = 0; index < requirement.size(); ++index) {
        if (coverage.sources[index] != ImpedanceSource::NOT_JUDGED) {
            judgedFrequencies.push_back(requirement[index].get_frequency());
        }
    }
    return judgedFrequencies;
}

std::optional<std::vector<double>> MagneticFilterImpedance::get_measured_frequencies(Magnetic* magnetic, Inputs* inputs) const {
    if (!inputs->get_design_requirements().get_minimum_impedance()) {
        return std::nullopt;
    }
    auto coverage = classify_requirement(magnetic, inputs);
    const auto requirement = inputs->get_design_requirements().get_minimum_impedance().value();
    std::vector<double> measuredFrequencies;
    for (size_t index = 0; index < requirement.size(); ++index) {
        if (coverage.sources[index] == ImpedanceSource::MEASURED) {
            measuredFrequencies.push_back(requirement[index].get_frequency());
        }
    }
    return measuredFrequencies;
}

std::pair<bool, double> MagneticFilterImpedance::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;

    // Impedance scoring: dimensionless log-ratio.
    //
    // For a "minimum impedance" requirement the part is invalid if Zact < Zreq at any
    // judged frequency. Among valid parts we want to reward parts close to the requirement
    // (smallest over-dimensioning) and penalize gross over-dimensioning (cost/size).
    //
    // Per-frequency:
    //     dev_i = log10(max(Zact_i / Zreq_i, 1))   // 0 when at-spec (under-spec → invalid)
    // Filter score = mean over the judged frequency points (then min-max normalized + inverted
    // downstream so smallest dev becomes the top score).
    //
    // No dead band on the raw score: normalize_scoring rescales to [0,1] so any
    // constant offset is wiped out, and a hard zero-clamp collapses all in-band parts
    // to the same value -- which then triggers the min==max degenerate branch in
    // normalize_scoring and returns 1.0 for every part (the "all show 100" bug).
    // Keeping the raw monotonic dev preserves ranking spread when this is the only
    // active filter.
    //
    // Judged frequencies: only the requirement points inside the core material's tabulated
    // complex-permeability range (get_judged_frequencies). Outside it the material has no mu(f),
    // so |Z| there is unknown: such a point is not judged at all -- neither a pass nor a fail,
    // and nothing is invented for it. A CISPR-band request reaches 30 MHz while many CMC ferrites
    // are tabulated only to 1-16 MHz; judging those parts on the band their data covers keeps
    // them in the ranking (El Choker: 142 of 299 chokes were excluded), and the mean over the
    // judged points keeps their score on the same per-point scale as a fully covered part's. The
    // adviser reports the judged frequencies per candidate. A part with no judged point at all
    // cannot be judged on the requirement and throws RequirementOutsideMaterialDataException.
    //
    // Measured band: a requirement point outside the material's mu(f) range is judged on the part's
    // own measured common-mode |Z| when its datasheet carries one (get_measured_impedance_curve:
    // the zero-bias commonModeChoke impedancePoints), interpolated log|Z| against log f between the
    // two bracketing measured points and never extrapolated past them. It is compared to the
    // requirement exactly as the model |Z| is (same validity, same per-point score). ACME's A07/A05
    // mu(f) ends at 12.5-14 MHz (their plot stops at mu = 10) while the WE CMCs on them are measured
    // to 1 GHz, so a 30 MHz CISPR point is judged on the measurement rather than left out. The
    // adviser reports these points as the candidate's measured frequencies.

    // Candidate SCORING runs the fast (OneLayer) capacitance path, explicitly. MKF d424c32e made
    // the full energy-based capacitance model Impedance's default -- right for analysing ONE
    // magnetic (that is the flagship path the Sweeper and the impedance panels use) -- but this
    // filter is called once per candidate for every core the adviser considers, and the full
    // model's per-turn energy sum costs seconds per candidate: the DMC default-wizard repro
    // (Test_CoreAdviser_DMC_Default_Wizard_Hang_Repro, 10 s budget) went to 28 minutes and the
    // INTERFERENCE_SUPPRESSION ranking snapshot moved 0.26 %. Ranking needs the same fast model
    // it was pinned and budgeted against; the chosen design is then analysed with the full one.
    constexpr bool kFastCapacitanceForScoring = true;

    if (inputs->get_design_requirements().get_minimum_impedance()) {
        auto impedanceRequirement = inputs->get_design_requirements().get_minimum_impedance().value();
        auto coverage = classify_requirement(magnetic, inputs);
        size_t numberJudged = std::count_if(coverage.sources.begin(), coverage.sources.end(),
                                            [](ImpedanceSource source) { return source != ImpedanceSource::NOT_JUDGED; });
        if (numberJudged == 0) {
            auto material = magnetic->get_core().resolve_material();
            auto [lowest, highest] = std::minmax_element(impedanceRequirement.begin(), impedanceRequirement.end(),
                [](const ImpedanceAtFrequency& a, const ImpedanceAtFrequency& b) { return a.get_frequency() < b.get_frequency(); });
            std::optional<std::pair<double, double>> measuredRange;
            if (coverage.measuredCurve) {
                measuredRange = std::make_pair(coverage.measuredCurve->front().first, coverage.measuredCurve->back().first);
            }
            throw RequirementOutsideMaterialDataException(material.get_name(), "minimumImpedance",
                                                          lowest->get_frequency(), highest->get_frequency(),
                                                          coverage.minimumMaterialFrequency, coverage.maximumMaterialFrequency,
                                                          measuredRange);
        }
        for (size_t index = 0; index < impedanceRequirement.size(); ++index) {
            const auto& impedanceAtFrequency = impedanceRequirement[index];
            if (coverage.sources[index] == ImpedanceSource::NOT_JUDGED) {
                continue;  // outside the material's mu(f) data and the measured curve: not judged
            }
            double zAct = coverage.sources[index] == ImpedanceSource::MODEL
                ? abs(OpenMagnetics::Impedance(kFastCapacitanceForScoring).calculate_impedance(*magnetic, impedanceAtFrequency.get_frequency()))
                : interpolate_measured_impedance(coverage.measuredCurve.value(), impedanceAtFrequency.get_frequency());
            double zReq = impedanceAtFrequency.get_impedance().get_magnitude();

            if (zReq > zAct) {
                valid = false;
            }

            // Ratio is clamped to >= 1: under-spec parts are already invalidated above,
            // so we only score over-dimensioning.
            double ratio = (zReq > 0) ? std::max(zAct / zReq, 1.0) : 1.0;
            double dev = std::log10(ratio);
            scoring += dev;
        }
        scoring /= numberJudged;
    }

    // Emit the impedance output at the operating points so downstream UI has the simulated |Z|.
    // We deliberately do NOT add this into `scoring`: there is no requirement here, and mixing
    // 1/|Z| in arbitrary units corrupts the min-max normalization (kills monotonicity and
    // frequency-fairness).
    //
    // This |Z| is display only, so it is written only at operating points whose frequency lies
    // inside the core material's tabulated complex-permeability span. Outside it there is no
    // mu(f) to evaluate (get_complex_permeability throws), and the output is left ABSENT at that
    // operating point rather than invented. A common-mode choke's inputs carry the 50 Hz mains
    // point (line current, thermal) next to the noise point, and every CMC ferrite is tabulated
    // from above 50 Hz: the throw here made the adviser drop the part over a number it only
    // shows (El Choker: 4 of 299 chokes left). The minimumImpedance requirement above is judged
    // inside the span only, and throws when no requirement point lies inside it.
    if (inputs->get_operating_points().size() > 0 && outputs != nullptr) {
        auto [minimumMaterialFrequency, maximumMaterialFrequency] = ComplexPermeability().get_frequency_range(magnetic->get_core().resolve_material());
        for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
            const auto& operatingPoint = inputs->get_operating_points()[operatingPointIndex];
            double frequency = operatingPoint.get_excitations_per_winding()[0].get_frequency();

            while (outputs->size() < operatingPointIndex + 1) {
                outputs->push_back(Outputs());
            }
            if (!(frequency >= minimumMaterialFrequency && frequency <= maximumMaterialFrequency)) {
                continue;
            }

            auto impedance = OpenMagnetics::Impedance(kFastCapacitanceForScoring).calculate_impedance(*magnetic, frequency);
            std::string name = magnetic->get_coil().get_functional_description()[0].get_name();
            ImpedanceOutput impedanceOutput;
            ComplexMatrixAtFrequency complexMatrixAtFrequency;
            complexMatrixAtFrequency.set_frequency(frequency);
            complexMatrixAtFrequency.get_mutable_magnitude()[name][name].set_nominal(abs(impedance));
            std::vector<ComplexMatrixAtFrequency> impedanceMatrixPerFrequency;
            impedanceMatrixPerFrequency.push_back(complexMatrixAtFrequency);
            impedanceOutput.set_impedance_matrix(impedanceMatrixPerFrequency);
            impedanceOutput.set_origin(ResultOrigin::SIMULATION);
            (*outputs)[operatingPointIndex].set_impedance(impedanceOutput);
        }
    }

    return {valid, scoring};
}

} // namespace OpenMagnetics
