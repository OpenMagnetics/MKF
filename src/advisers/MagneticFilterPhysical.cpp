#include "advisers/MagneticFilter.h"
#include "advisers/MagneticFilterInternal.h"
#include "physical_models/MagnetizingInductance.h"
#include "processors/Inputs.h"
#include "processors/ResultsAnalyzer.h"
#include "support/Exceptions.h"
#include "support/Settings.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace OpenMagnetics {

// Phase 5: extracted from MagneticFilter.cpp.
// This translation unit owns the saturation / current-density filters:
//   - MagneticFilterSaturation (flux density and saturation current from the magnetizing current)
//   - MagneticFilterDcCurrentDensity
//   - MagneticFilterEffectiveCurrentDensity
// The shared magnetizing-current helpers live in MagneticFilterInternal.h.

std::pair<bool, double> MagneticFilterSaturation::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;

    // One saturation check for every design, without classifying it: the flux density is the one
    // the magnetizing-inductance model gives at the operating point, from the magnetizing current it
    // stamps (on a single winding the winding current; with several windings the magnetizing current
    // of the circuit: Faraday on the voltage, which is the volt-seconds' flux whatever the
    // permeability, or the ampere-turns of the currents). The saturation-current gate compares the
    // gap-aware saturation current with that same magnetizing current, so a transformer is never
    // gated on its load current and a coupled inductor never escapes the gate.
    const bool singleWinding = is_single_winding_design(*inputs);

    const std::string magneticRef = inductance_flux_cache_key(*magnetic);
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
        // Evaluate B_sat at the hot junction corner (ABT #13), not the 25 C
        // ambient: ferrite B_sat falls with temperature, so the hot corner is the
        // real operating point. RAW B_sat (proportion=false) — the multiplicative
        // `margin` below is the only safety factor on this gate; the 0.7
        // flux-proportion is a SIZING target, not stacked here.
        auto magneticFluxDensitySaturation = magnetic->get_mutable_core().get_magnetic_flux_density_saturation(
            saturation_derating_temperature(operatingPoint.get_conditions().get_ambient_temperature()),
            /*proportion=*/false);
        
        // The DC point of the magnetizing model: the saturation check's second criterion
        // (the bias collapsed the permeability to mu_r = 1, Alf 2026-10-14).
        std::optional<DcBiasPoint> dcBiasPoint;
        {
            OpenMagnetics::MagnetizingInductance magnetizingInductanceObj;

            magnetizingInductanceObj.set_topology(magnetizing_model_topology(*inputs));
            // The model gives its numbers for a saturated core too (ABT #1855): the verdict is the saturation
            // check below on the flux density it returns, as analyze_results makes it.
            std::pair<MagnetizingInductanceOutput, SignalDescriptor> result =
                magnetizingInductanceObj.calculate_inductance_and_magnetic_flux_density(
                    magnetic->get_core(), magnetic->get_coil(), &operatingPoint);
            auto magneticFluxDensity = result.second;
            dcBiasPoint = magnetizingInductanceObj.get_last_dc_bias_point();

            // Phase 8 (perf): seed the per-adviser-invocation cache so
            // downstream filters (TEMPERATURE) that need the same B for
            // the same (magnetic, OP) can skip the iterative recomputation.
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

        // Common-mode chokes have every winding on one isolation side, but their line current is DIFFERENTIAL: the flux of
        // the line/neutral windings cancels in the core, so the raw line
        // current does NOT saturate it. The only flux driver is the CM noise
        // current, and that is NOT DERIVABLE from the design inputs: a
        // symmetric operating point stores the line/ripple current identically
        // on every winding (DM, cancels), and MAS has no CM-noise-current
        // requirement to gate against. get_common_mode_choke_magnetizing_current
        // models the ripple swing AS the CM current — an over-estimate by
        // construction (the ripple returns on the other conductor and cancels
        // too), which put multi-tesla fictional B on every EMI toroid and
        // zeroed the default web CMC design (ABT #551/#236 archaeology: the
        // gate rejected all 25 impedance survivors, so the test never passed
        // with its committed fixture). BOTH saturation gates therefore skip
        // tagged CMCs: saturation for a CMC is only checkable against an
        // explicit CM-current spec, which does not exist in MAS today.
        //
        // CMC-ONLY, deliberately: a differential-mode choke ALSO has 2-4
        // matching-current windings (can_be_common_mode_choke is true for it
        // too), but a DMC's flux does NOT cancel — it carries the full line
        // current additively and genuinely saturates — so a DMC must KEEP both
        // gates. Distinguish by the CMC tag (topology or the web wizard's
        // subApplication), then confirm current symmetry with
        // can_be_common_mode_choke so a mis-tagged asymmetric design still
        // gets checked.
        // Shared with the MagneticAdviser's final gate (MagneticFilterInternal.h).
        bool isCommonModeChoke = is_tagged_common_mode_choke(*inputs, operatingPoint);

        // Saturation scoring: dimensionless ratio of operating peak flux density to
        // material saturation. Smaller is better (more headroom). The previous
        // implementation used fabs(Bsat - Bpeak), which under min-max+invert
        // *rewards* candidates that operate close to saturation - the opposite of
        // engineering intent. Ratio Bpeak/Bsat preserves the "more headroom = higher
        // normalized score" semantics correctly. For a CMC the B here is the
        // fictional one described above, so neither score nor gate uses it.
        const double margin = Settings::GetInstance().get_core_adviser_saturation_margin();
        if (!isCommonModeChoke) {
            double bRatio = (magneticFluxDensitySaturation > 0)
                ? (magneticFluxDensityPeak / magneticFluxDensitySaturation)
                : 1.0;
            scoring += bRatio;

            // Multiplicative derating: reject if Bpeak·margin > Bsat (i.e. the
            // core would saturate at `margin`× the actual operating excursion).
            // margin == 1.0 (default) reproduces the historical bare comparison.
            // Set via Settings::set_core_adviser_saturation_margin() — Maniktala
            // Ch.5 recommends ≥1.2 for ferrite designs to leave headroom for
            // tolerance, temperature, and DC-bias swing. The check is the one
            // analyze_results makes (ResultsAnalyzer::check_saturation, ABT #1855),
            // here with the adviser's margin and derating temperature.
            auto saturation = ResultsAnalyzer::check_saturation(
                magnetic->get_mutable_core().resolve_material(), magneticFluxDensityPeak,
                saturation_derating_temperature(operatingPoint.get_conditions().get_ambient_temperature()), margin, dcBiasPoint);

            if (saturation.saturates) {
                logEntry(std::string("Saturation filter: rejecting candidate: ") + saturation.message, "MagneticFilter", 2);
                return {false, 0.0};
            }
        }

        if (!isCommonModeChoke) {
            // Authoritative saturation-current gate: require the gap-aware
            // saturation current to clear the peak magnetizing current by the
            // margin. The B-based check above evaluates B at the operating current
            // with permeability rolloff, which is LOOSER for gapped inductors:
            // candidates passing it could still saturate (canonical buck
            // 36/48/60 V→12 V/5 A/200 kHz: EQ 13/6 passed on B but had isat
            // 2.93 A < 1.2 × ipeak 3.19 A).
            //
            // The current is the magnetizing current: on a single winding its own
            // current; with several windings the one the model stamped above (for a
            // transformer the volt-seconds' magnetizing current, never the load
            // current the windings exchange; for a coupled inductor or flyback the
            // net ampere-turns or the Faraday current of the energy it stores).
            //
            // Derating (ABT #13): evaluate I_sat at the hot junction corner and
            // with RAW B_sat (proportion=false). The total derating on this gate
            // is exactly (hot temperature × margin) — the 0.7 flux-proportion that
            // calculate_saturation_current applies by default is a SIZING target
            // and is deliberately NOT stacked here, so temperature + margin don't
            // compound into a runaway combined margin.
            //
            // Pool starvation (the failure mode of the reverted hard reject
            // 2181dd1c) is handled by the pipeline retry, which relaxes the
            // margin to 1.0 and tags the results.
            double currentPeak;
            if (singleWinding) {
                auto excitation = Inputs::get_primary_excitation(operatingPoint);
                if (!excitation.get_current() || !excitation.get_current()->get_processed() ||
                    !excitation.get_current()->get_processed()->get_peak()) {
                    throw InvalidInputException(ErrorCode::MISSING_DATA,
                        "Saturation filter needs the processed peak of the primary current to gate by saturation current");
                }
                currentPeak = excitation.get_current()->get_processed()->get_peak().value();
            }
            else {
                currentPeak = stamped_magnetizing_current_peak(operatingPoint);
            }
            double saturationCurrent = magnetic->calculate_saturation_current(
                saturation_derating_temperature(operatingPoint.get_conditions().get_ambient_temperature()),
                /*proportion=*/false);

            // Headroom utilization (margin·ipeak/isat): < 1 clears the gate,
            // smaller = more headroom — same direction as bRatio above.
            scoring += (saturationCurrent > 0) ? (margin * currentPeak / saturationCurrent) : 1.0;

            if (saturationCurrent < margin * currentPeak) {
                return {false, 0.0};
            }
        }

        ++opIndex;
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
