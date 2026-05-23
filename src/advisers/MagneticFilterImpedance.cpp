#include "advisers/MagneticFilter.h"
#include "advisers/MagneticFilterInternal.h"
#include "constructive_models/NumberTurns.h"
#include "physical_models/Impedance.h"
#include "physical_models/WindingLosses.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "support/Exceptions.h"

#include <algorithm>
#include <cmath>
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
                    zUnit = abs(_impedanceModel.calculate_impedance(core, unitCoil, frequency, H_dc_unit, defaults.ambientTemperature));
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
                            impedance = abs(_impedanceModel.calculate_impedance(core, coil, frequency, H_dc, defaults.ambientTemperature));
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

std::pair<bool, double> MagneticFilterImpedance::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;

    // Impedance scoring: dimensionless log-ratio.
    //
    // For a "minimum impedance" requirement the part is invalid if Zact < Zreq at any
    // frequency. Among valid parts we want to reward parts close to the requirement
    // (smallest over-dimensioning) and penalize gross over-dimensioning (cost/size).
    //
    // Per-frequency:
    //     dev_i = log10(max(Zact_i / Zreq_i, 1))   // 0 when at-spec (under-spec → invalid)
    // Filter score = mean over frequency points (then min-max normalized + inverted
    // downstream so smallest dev becomes the top score).
    //
    // No dead band on the raw score: normalize_scoring rescales to [0,1] so any
    // constant offset is wiped out, and a hard zero-clamp collapses all in-band parts
    // to the same value -- which then triggers the min==max degenerate branch in
    // normalize_scoring and returns 1.0 for every part (the "all show 100" bug).
    // Keeping the raw monotonic dev preserves ranking spread when this is the only
    // active filter.

    if (inputs->get_design_requirements().get_minimum_impedance()) {
        auto impedanceRequirement = inputs->get_design_requirements().get_minimum_impedance().value();
        for (auto impedanceAtFrequency : impedanceRequirement) {
            auto impedance = OpenMagnetics::Impedance().calculate_impedance(*magnetic, impedanceAtFrequency.get_frequency());
            double zReq = impedanceAtFrequency.get_impedance().get_magnitude();
            double zAct = abs(impedance);

            if (zReq > zAct) {
                valid = false;
            }

            // Ratio is clamped to >= 1: under-spec parts are already invalidated above,
            // so we only score over-dimensioning.
            double ratio = (zReq > 0) ? std::max(zAct / zReq, 1.0) : 1.0;
            double dev = std::log10(ratio);
            scoring += dev;
        }
        scoring /= impedanceRequirement.size();
    }

    // Always emit the impedance output for any operating points so downstream UI
    // has the simulated |Z|. We deliberately do NOT add this into `scoring`: there
    // is no requirement here, and mixing 1/|Z| in arbitrary units corrupts the
    // min-max normalization (kills monotonicity and frequency-fairness).
    if (inputs->get_operating_points().size() > 0 && outputs != nullptr) {
        for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
            auto operatingPoint = inputs->get_operating_points()[operatingPointIndex];
            auto impedance = OpenMagnetics::Impedance().calculate_impedance(*magnetic, operatingPoint.get_excitations_per_winding()[0].get_frequency());
            std::string name = magnetic->get_coil().get_functional_description()[0].get_name();

            while (outputs->size() < operatingPointIndex + 1) {
                outputs->push_back(Outputs());
            }
            ImpedanceOutput impedanceOutput;
            ComplexMatrixAtFrequency complexMatrixAtFrequency;
            complexMatrixAtFrequency.set_frequency(operatingPoint.get_excitations_per_winding()[0].get_frequency());
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
