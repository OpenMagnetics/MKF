#include "advisers/MagneticFilter.h"
#include "advisers/MagneticFilterInternal.h"
#include "physical_models/MagnetizingInductance.h"
#include "processors/Inputs.h"
#include "support/Exceptions.h"
#include "support/Settings.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace OpenMagnetics {

// Phase 5: extracted from MagneticFilter.cpp.
// This translation unit owns the saturation / current-density filters:
//   - MagneticFilterSaturation (transformer-vs-inductor B detection)
//   - MagneticFilterDcCurrentDensity
//   - MagneticFilterEffectiveCurrentDensity
// The is_energy_storing_topology() helper lives in MagneticFilterInternal.h
// and is shared with the (still-monolithic) MagneticFilter.cpp.

std::pair<bool, double> MagneticFilterSaturation::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;

    // Transformer vs Inductor Detection
    // ==================================
    // Transformers: B is determined by applied VOLTAGE (Faraday's Law)
    //   - B_peak = V_peak / (N × Ae × ω)
    //   - No permeability iteration needed since magnetizing current is ideally small
    //
    // Inductors/Energy-storing: B is determined by CURRENT and permeability
    //   - Uses iterative calculation accounting for permeability rolloff with DC bias
    //
    // Detection priority:
    // 1. Use topology if specified (most reliable)
    // 2. Fall back to inductance field heuristic (minimum-only = transformer)
    //
    auto topology = inputs->get_design_requirements().get_topology();
    bool isTransformer;
    if (windings_on_single_isolation_side(inputs->get_design_requirements().get_isolation_sides())) {
        // All windings on one isolation side -> (coupled) inductor, never a transformer,
        // regardless of the converter topology (e.g. Weinberg L1 input coupled inductor).
        // Forces the B-from-current (DC-biased) path instead of B-from-voltage.
        isTransformer = false;
    } else if (topology.has_value()) {
        // Use topology-based detection
        isTransformer = !is_energy_storing_topology(topology);
    } else {
        // Legacy heuristic: minimum-only inductance = transformer
        isTransformer = inputs->get_design_requirements().get_magnetizing_inductance().get_minimum() &&
                         !inputs->get_design_requirements().get_magnetizing_inductance().get_nominal() &&
                         !inputs->get_design_requirements().get_magnetizing_inductance().get_maximum();
    }

    const std::string magneticRef = magnetic->get_reference();
    size_t opIndex = 0;
    for (auto operatingPoint : inputs->get_operating_points()) {
        double magneticFluxDensityPeak;
        // Use RAW B_sat here (proportion=false). The Maniktala-style
        // multiplicative `margin` derating below (Settings::get_core_adviser_-
        // saturation_margin(), default 1.2) is the saturation-side safety
        // factor. The `proportion` factor — defaults.maximumProportionMagnetic-
        // FluxDensitySaturation = 0.7 — is the design-side target headroom
        // used at initial-turn selection time (CoreAdviserDataset.cpp:551)
        // so the picked N keeps steady-state B_peak ≈ 70 % of B_sat_raw.
        // Stacking both factors here was a double-count: the filter rejected
        // any candidate where B_peak × 1.2 > 0.7 × B_sat_raw (i.e.
        // B_peak > 0.58 × B_sat_raw), which the initial-N seed cannot satisfy
        // because it deliberately picked N for B_peak ≈ 0.7 × B_sat_raw.
        // Net effect was every PSFB / LLC FB / SRC transformer that depended
        // on the saturation-aware seed losing its entire candidate set at
        // this filter (HANDOFF_HEAVY_TEST_GAPS.md §Gap 3 PSFB diagnosis).
        auto magneticFluxDensitySaturation = magnetic->get_mutable_core().get_magnetic_flux_density_saturation(
            operatingPoint.get_conditions().get_ambient_temperature(),
            /*proportion=*/false);
        
        if (isTransformer) {
            // For transformers, calculate B from voltage using MagnetizingInductance method
            // This avoids the permeability iteration issue with current-based calculation
            auto excitation = Inputs::get_primary_excitation(operatingPoint);
            double frequency = excitation.get_frequency();
            double voltagePeak = 0;
            double actualInductance = 0;
            
            bool hasVoltage = excitation.get_voltage() && excitation.get_voltage()->get_processed() && 
                excitation.get_voltage()->get_processed()->get_peak();
            bool hasCurrentWaveform = excitation.get_current() && excitation.get_current()->get_waveform();
            
            if (hasVoltage) {
                voltagePeak = excitation.get_voltage()->get_processed()->get_peak().value();
            } else if (hasCurrentWaveform) {
                // Derive voltage from V = L·di/dt using the actual gapped-core inductance
                OpenMagnetics::MagnetizingInductance magnetizingInductanceCalc;
                auto inductanceOutput = magnetizingInductanceCalc.calculate_inductance_from_number_turns_and_gapping(
                    *magnetic, &operatingPoint);

                if (inductanceOutput.get_magnetizing_inductance().get_nominal()) {
                    actualInductance = inductanceOutput.get_magnetizing_inductance().get_nominal().value();
                    if (actualInductance > 0) {
                        auto inducedVoltage = Inputs::calculate_induced_voltage(excitation, actualInductance);
                        if (inducedVoltage.get_processed() && inducedVoltage.get_processed()->get_peak()) {
                            voltagePeak = inducedVoltage.get_processed()->get_peak().value();
                        }
                    }
                }
            }
            
            double numberTurns = magnetic->get_coil().get_functional_description()[0].get_number_turns();
            OpenMagnetics::MagnetizingInductance magnetizingInductanceObj;

            // For non-sinusoidal waveforms (DAB square, CLLLC trapezoid, ...)
            // the sinusoidal V_peak/(N·A_e·ω) underestimates B_peak by ~36 %.
            // Use the V·s excursion path (Faraday integrated).
            double maxVoltSeconds = hasVoltage ? Inputs::calculate_max_volt_seconds(excitation) : 0.0;
            if (maxVoltSeconds > 0 && numberTurns > 0) {
                magneticFluxDensityPeak = magnetizingInductanceObj
                    .calculate_flux_density_peak_from_volt_seconds(
                        magnetic->get_mutable_core(), numberTurns, maxVoltSeconds);
            } else {
                magneticFluxDensityPeak = magnetizingInductanceObj.calculate_flux_density_peak_from_voltage(
                    magnetic->get_mutable_core(), numberTurns, voltagePeak, frequency);
            }
        } else {
            // For inductors/energy-storing converters, use the standard calculation
            OpenMagnetics::MagnetizingInductance magnetizingInductanceObj;
            
            auto result = magnetizingInductanceObj.calculate_inductance_and_magnetic_flux_density(
                magnetic->get_core(), magnetic->get_coil(), &operatingPoint);
            auto magneticFluxDensity = result.second;

            // Phase 8 (perf): seed the per-adviser-invocation cache so
            // downstream filters (TEMPERATURE) that need the same B for
            // the same (magnetic, OP) can skip the iterative recomputation.
            // Inductor path only — transformer branch uses a different B
            // derivation and shouldn't share the cache.
            if (magneticFluxDensity.get_processed() && magneticFluxDensity.get_processed()->get_peak()
                && result.first.get_magnetizing_inductance().get_nominal()) {
                cache_inductance_flux(magneticRef, opIndex,
                    result.first.get_magnetizing_inductance().get_nominal().value(),
                    magneticFluxDensity);
            }
            
            if (!magneticFluxDensity.get_processed() || !magneticFluxDensity.get_processed()->get_peak()) {
                return {false, 0.0};
            }
            magneticFluxDensityPeak = magneticFluxDensity.get_processed()->get_peak().value();
        }

        // Saturation scoring: dimensionless ratio of operating peak flux density to
        // material saturation. Smaller is better (more headroom). The previous
        // implementation used fabs(Bsat - Bpeak), which under min-max+invert
        // *rewards* candidates that operate close to saturation - the opposite of
        // engineering intent. Ratio Bpeak/Bsat preserves the "more headroom = higher
        // normalized score" semantics correctly.
        double bRatio = (magneticFluxDensitySaturation > 0)
            ? (magneticFluxDensityPeak / magneticFluxDensitySaturation)
            : 1.0;
        scoring += bRatio;

        // Multiplicative derating: reject if Bpeak·margin > Bsat (i.e. the
        // core would saturate at `margin`× the actual operating excursion).
        // margin == 1.0 (default) reproduces the historical bare comparison.
        // Set via Settings::set_core_adviser_saturation_margin() — Maniktala
        // Ch.5 recommends ≥1.2 for ferrite designs to leave headroom for
        // tolerance, temperature, and DC-bias swing.
        const double margin = Settings::GetInstance().get_core_adviser_saturation_margin();
        bool isSaturated = magneticFluxDensityPeak * margin > magneticFluxDensitySaturation;

        ++opIndex;

        if (isSaturated) {
            return {false, 0.0};
        }
    }

    scoring /= inputs->get_operating_points().size();

    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterDcCurrentDensity::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;

    if (inputs->get_operating_points().size() > 0 && magnetic->get_mutable_coil().get_functional_description().size() != inputs->get_operating_points()[0].get_excitations_per_winding().size()) {
        return {false, 0.0};        
    }

    for (const auto& operatingPoint : inputs->get_operating_points()) {
        for (size_t windingIndex = 0; windingIndex < magnetic->get_mutable_coil().get_functional_description().size(); ++windingIndex) {
            auto excitation = operatingPoint.get_excitations_per_winding()[windingIndex];
            if (!excitation.get_current()) { 
                throw InvalidInputException(ErrorCode::MISSING_DATA, "Current is missing in excitation");
            }
            auto current = excitation.get_current().value();
            auto wire = magnetic->get_mutable_coil().resolve_wire(windingIndex);
            auto dcCurrentDensity = wire.calculate_dc_current_density(current) / magnetic->get_mutable_coil().get_functional_description()[windingIndex].get_number_parallels();

            // DC current density scoring: ratio of actual to maximum allowed.
            // Smaller is better (more headroom against thermal limit). See note on
            // MagneticFilterSaturation for why fabs(max - actual) was wrong.
            scoring += dcCurrentDensity / defaults.maximumCurrentDensity;
            if (dcCurrentDensity > defaults.maximumCurrentDensity) {
                return {false, 0.0};
            }
        }
    }

    scoring /= inputs->get_operating_points().size();

    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterEffectiveCurrentDensity::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;

    if (inputs->get_operating_points().size() > 0 && magnetic->get_mutable_coil().get_functional_description().size() != inputs->get_operating_points()[0].get_excitations_per_winding().size()) {
        return {false, 0.0};        
    }

    for (const auto& operatingPoint : inputs->get_operating_points()) {
        for (size_t windingIndex = 0; windingIndex < magnetic->get_mutable_coil().get_functional_description().size(); ++windingIndex) {
            auto excitation = operatingPoint.get_excitations_per_winding()[windingIndex];
            if (!excitation.get_current()) {
                throw InvalidInputException(ErrorCode::MISSING_DATA, "Current is missing in excitation");
            }
            auto current = excitation.get_current().value();
            auto wire = magnetic->get_mutable_coil().resolve_wire(windingIndex);
            auto effectiveCurrentDensity = wire.calculate_effective_current_density(current, operatingPoint.get_conditions().get_ambient_temperature()) / magnetic->get_mutable_coil().get_functional_description()[windingIndex].get_number_parallels();

            // Effective (AC+DC, temperature-corrected) current density scoring:
            // ratio of actual to maximum allowed. Smaller is better. See note on
            // MagneticFilterSaturation for the rationale.
            scoring += effectiveCurrentDensity / defaults.maximumEffectiveCurrentDensity;
            if (effectiveCurrentDensity > defaults.maximumEffectiveCurrentDensity) {
                return {false, 0.0};
            }
        }
    }

    scoring /= inputs->get_operating_points().size();

    return {valid, scoring};
}

} // namespace OpenMagnetics
