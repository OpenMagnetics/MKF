// Phase 5: loss-filter implementations extracted from MagneticFilter.cpp.
// All four classes (CoreAndDcLosses, CoreDcAndSkinLosses, Losses,
// LossesNoProximity) share the helpers in MagneticFilterInternal.h and
// nothing else, so they sit naturally in their own translation unit.
//
// Class declarations remain in advisers/MagneticFilter.h; this file
// supplies the definitions only.

#include "advisers/MagneticFilter.h"
#include "advisers/MagneticFilterInternal.h"
#include "constructive_models/NumberTurns.h"
#include "physical_models/WindingLosses.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "support/Exceptions.h"
#include "support/Utils.h"

#include <cmath>
#include <limits>

namespace OpenMagnetics {

MagneticFilterCoreAndDcLosses::MagneticFilterCoreAndDcLosses(Inputs inputs)
    : MagneticFilterCoreAndDcLosses(inputs, default_loss_filter_models()) {}

MagneticFilterCoreAndDcLosses::MagneticFilterCoreAndDcLosses() {
    auto models = default_loss_filter_models();
    _maximumPowerMean = 0;
    _coreLossesModelSteinmetz = CoreLossesModel::factory(models);
    _coreLossesModelProprietary = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Proprietary"}}));
    _magnetizingInductance = MagnetizingInductance(models["gapReluctance"]);
    _models = models;
}

MagneticFilterCoreAndDcLosses::MagneticFilterCoreAndDcLosses(Inputs inputs, std::map<std::string, std::string> models) {
    _maximumPowerMean = compute_maximum_power_mean_and_maybe_force_steinmetz(inputs, models);
    _coreLossesModelSteinmetz = CoreLossesModel::factory(models);
    _coreLossesModelProprietary = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Proprietary"}}));
    _magnetizingInductance = MagnetizingInductance(models["gapReluctance"]);
    _models = models;
}

std::pair<bool, double> MagneticFilterCoreAndDcLosses::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    const auto& core = magnetic->get_core();

    if (inputs->get_operating_points().size() > 0 && magnetic->get_mutable_coil().get_functional_description().size() != inputs->get_operating_points()[0].get_excitations_per_winding().size()) {
        return {false, 0.0};
    }

    // This filter winds candidates with fast_wind() (no delimit/compact) for
    // speed. The flag lives on the process-global Settings singleton, so it must
    // be restored on every exit path — including the many throwing/early-return
    // paths below — or it leaks into all later code that reads it. RAII instead
    // of manual save/restore (which previously had no restore at all).
    SettingsGuard<bool> coilDelimitGuard(settings,
        &Settings::get_coil_delimit_and_compact,
        &Settings::set_coil_delimit_and_compact, false);

    std::string shapeName = core.get_shape_name();
    prepare_bobbin_for_non_pqi(magnetic, shapeName);

    auto currentNumberTurns = magnetic->get_coil().get_functional_description()[0].get_number_turns();
    NumberTurns numberTurns(currentNumberTurns);
    std::vector<double> totalLossesPerOperatingPoint;
    std::vector<CoreLossesOutput> coreLossesPerOperatingPoint;
    std::vector<WindingLossesOutput> windingLossesPerOperatingPoint;
    double currentTotalLosses = std::numeric_limits<double>::infinity();
    double coreLosses = std::numeric_limits<double>::infinity();
    CoreLossesOutput coreLossesOutput;
    double ohmicLosses = std::numeric_limits<double>::infinity();
    WindingLossesOutput windingLossesOutput;
    windingLossesOutput.set_origin(ResultOrigin::SIMULATION);
    double newTotalLosses = std::numeric_limits<double>::infinity();
    auto previousNumberTurnsPrimary = currentNumberTurns;

    size_t iteration = defaults.coreAdviserSkinEffectMaxIterations;

    Coil coil = magnetic->get_coil();

    // Energy-storing inductors saturate as turns rise: at a fixed gap isat ∝ 1/N
    // (L ∝ N², isat = B_sat·N·A_e/L). The loss-minimising turn sweep below would
    // otherwise push N up (lower core loss) past the saturation-safe point — the
    // CoreAdviser's saturation-safe turn/gap sizing then gets silently undone and
    // the design saturates. Cap the sweep using the SAME
    // Magnetic::calculate_saturation_current the saturation filter and realism
    // gate use, so the loss-optimal turn count we return still clears the margin.
    bool isInductor = is_inductor(*inputs);
    double saturationMargin = Settings::GetInstance().get_core_adviser_saturation_margin();

    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
        auto operatingPoint = inputs->get_operating_point(operatingPointIndex);
        // Derating (ABT #13): hot junction corner, consistent with the saturation
        // filter's inductor gate (raw B_sat below).
        double temperature = saturation_derating_temperature(operatingPoint.get_conditions().get_ambient_temperature());
        OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];
        double saturationPeakCurrent = (isInductor && excitation.get_current() && excitation.get_current()->get_processed()
                                        && excitation.get_current()->get_processed()->get_peak())
                                       ? std::abs(excitation.get_current()->get_processed()->get_peak().value()) : 0.0;
        size_t numberTimeouts = 0;

        // Track the loss-minimum across the sweep. The do-while only gates progress;
        // it exits once losses STOP improving, so the coil state at break time is the
        // first WORSE point, not the optimum. Remember the best and restore it after
        // the loop (ABT #105 — ported from the CoreDcAndSkin twin, which already did
        // this; the CoreAndDc twin was pushing the last, worse iteration).
        double bestTotalLosses = std::numeric_limits<double>::infinity();
        uint64_t bestNumberTurnsPrimary = currentNumberTurns;
        CoreLossesOutput bestCoreLossesOutput;
        WindingLossesOutput bestWindingLossesOutput;
        bestWindingLossesOutput.set_origin(ResultOrigin::SIMULATION);

        do {
            currentTotalLosses = newTotalLosses;
            auto numberTurnsCombination = numberTurns.get_next_number_turns_combination();
            coil.get_mutable_functional_description()[0].set_number_turns(numberTurnsCombination[0]);
            // coil = Coil(coil);  delimit/compact disabled by coilDelimitGuard above
            coil.fast_wind();

            // Saturation cap for energy-storing inductors: this (higher) turn count
            // is checked against the gap-aware saturation current; once it dips
            // below margin·I_pk, every higher N saturates too, so revert to the last
            // safe turn count and stop sweeping.
            if (saturationPeakCurrent > 0) {
                Magnetic saturationMagnetic = *magnetic;
                saturationMagnetic.set_coil(coil);
                double saturationCurrent = 0;
                try { saturationCurrent = saturationMagnetic.calculate_saturation_current(temperature, /*proportion=*/false); }
                catch (const std::exception& e) {
                    // ABT #121.4: was silent — with saturationCurrent = 0 the
                    // saturation cap on the turns sweep is disabled, so say so.
                    logEntry(std::string("Losses filter: saturation current failed during turns sweep: ") +
                             e.what() + " — saturation cap disabled for this sweep step", "CoreAdviser", 2);
                    saturationCurrent = 0;
                }
                if (saturationCurrent > 0 && saturationCurrent < saturationMargin * saturationPeakCurrent) {
                    coil.get_mutable_functional_description()[0].set_number_turns(previousNumberTurnsPrimary);
                    settings.set_coil_delimit_and_compact(false);
                    coil.fast_wind();
                    break;
                }
            }

            auto [magnetizingInductance, magneticFluxDensity] = _magnetizingInductance.calculate_inductance_and_magnetic_flux_density(core, coil, &operatingPoint);

            if (!check_requirement(inputs->get_design_requirements().get_magnetizing_inductance(), magnetizingInductance.get_magnetizing_inductance().get_nominal().value())) {
                if (resolve_dimensional_values(inputs->get_design_requirements().get_magnetizing_inductance()) < resolve_dimensional_values(magnetizingInductance.get_magnetizing_inductance().get_nominal().value())) {
                    coil.get_mutable_functional_description()[0].set_number_turns(previousNumberTurnsPrimary);
                    // coil = Coil(coil);  delimit/compact disabled by coilDelimitGuard above
                    coil.fast_wind();
                    break;
                }
            }
            else {
                previousNumberTurnsPrimary = numberTurnsCombination[0];
            }

            if (!is_pqi_or_ui_shape(shapeName)) {
                if (!coil.get_turns_description()) {
                    // Phase 1 fix: previously this silently set
                    // newTotalLosses = coreLosses (still DBL_MAX on the
                    // first iteration) and broke out, relying on the
                    // OP-count mismatch later to fail the magnetic. Make
                    // the rejection explicit: fast_wind() could not lay
                    // out the coil, so this candidate is genuinely
                    // unusable (not a "core losses only" fallback).
                    return {false, 0.0};
                }
            }

            excitation.set_magnetic_flux_density(magneticFluxDensity);
            {
                auto pick = compute_core_losses_with_negative_guard(core, excitation, temperature, _coreLossesModelSteinmetz, _coreLossesModelProprietary);
                if (!pick.ok) {
                    // Phase 1 fix: was a silent `break` that left coreLosses
                    // negative and relied on a downstream guard. Proprietary
                    // models occasionally return tiny negatives (~µW) from
                    // interpolation noise at low excitation.
                    return {false, 0.0};
                }
                coreLossesOutput = pick.output;
                coreLosses = pick.value;
            }

            if (coreLosses < 0) {
                throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happend in core losses calculation for magnetic: " + magnetic->get_manufacturer_info().value().get_reference().value());
            }

            if (!coil.get_turns_description()) {
                // Phase 1 fix: same rationale as above — explicit reject
                // rather than the silent break + DBL_MAX leak.
                return {false, 0.0};
            }

            if (!is_pqi_or_ui_shape(shapeName)) {
                windingLossesOutput = _windingOhmicLosses.calculate_ohmic_losses(coil, operatingPoint, temperature);
                ohmicLosses = windingLossesOutput.get_winding_losses();
                newTotalLosses = coreLosses + ohmicLosses;
                if (ohmicLosses < 0) {
                    throw CalculationException(ErrorCode::CALCULATION_INVALID_INPUT, "Something wrong happend in ohmic losses calculation for magnetic: " + magnetic->get_manufacturer_info().value().get_reference().value() + " ohmicLosses: " + std::to_string(ohmicLosses));
                }
            }
            else {
                // PQI / UI shapes: integrated-winding layouts that fast_wind()
                // doesn't materialise, so ohmic losses can't be computed here.
                // Use core-only losses and stop sweeping turns. Explicit
                // per-shape policy — not a silent fallback.
                newTotalLosses = coreLosses;
                break;
            }

            if (!std::isfinite(newTotalLosses)) {
                throw CalculationException(ErrorCode::CALCULATION_DIVERGED, "Too large losses");
            }

            // Record the best point seen so far in this operating-point sweep.
            if (newTotalLosses < bestTotalLosses) {
                bestTotalLosses = newTotalLosses;
                bestNumberTurnsPrimary = numberTurnsCombination[0];
                bestCoreLossesOutput = coreLossesOutput;
                bestWindingLossesOutput = windingLossesOutput;
            }

            iteration--;
            if (iteration <=0) {
                numberTimeouts++;
                break;
            }
        }
        while(newTotalLosses < currentTotalLosses * defaults.coreAdviserThresholdValidity);

        // Restore the best N from the sweep so downstream code sees the loss-optimal coil.
        if (std::isfinite(bestTotalLosses) &&
            coil.get_functional_description()[0].get_number_turns() != static_cast<double>(bestNumberTurnsPrimary)) {
            coil.get_mutable_functional_description()[0].set_number_turns(bestNumberTurnsPrimary);
            coil.fast_wind();  // delimit/compact disabled by coilDelimitGuard above
        }

        if (std::isfinite(bestTotalLosses)) {
            magnetic->set_coil(coil);
            totalLossesPerOperatingPoint.push_back(bestTotalLosses);
            coreLossesPerOperatingPoint.push_back(bestCoreLossesOutput);
            windingLossesPerOperatingPoint.push_back(bestWindingLossesOutput);
        }
        else if (std::isfinite(coreLosses) && coreLosses > 0) {
            // Fallback: no point ever became finite (e.g. PQI/UI shortcut path with only core losses).
            magnetic->set_coil(coil);
            currentTotalLosses = newTotalLosses;
            totalLossesPerOperatingPoint.push_back(currentTotalLosses);
            coreLossesPerOperatingPoint.push_back(coreLossesOutput);
            windingLossesPerOperatingPoint.push_back(windingLossesOutput);
        }
    }

    return finalize_losses_scoring(totalLossesPerOperatingPoint, coreLossesPerOperatingPoint, windingLossesPerOperatingPoint,
                                    magnetic, inputs, outputs, _maximumPowerMean);
}

MagneticFilterCoreDcAndSkinLosses::MagneticFilterCoreDcAndSkinLosses(Inputs inputs)
    : MagneticFilterCoreDcAndSkinLosses(inputs, default_loss_filter_models()) {}

MagneticFilterCoreDcAndSkinLosses::MagneticFilterCoreDcAndSkinLosses() {
    auto models = default_loss_filter_models();
    _maximumPowerMean = 0;
    _coreLossesModelSteinmetz = CoreLossesModel::factory(models);
    _coreLossesModelProprietary = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Proprietary"}}));
    _magnetizingInductance = MagnetizingInductance(models["gapReluctance"]);
    _models = models;
}

MagneticFilterCoreDcAndSkinLosses::MagneticFilterCoreDcAndSkinLosses(Inputs inputs, std::map<std::string, std::string> models) {
    _maximumPowerMean = compute_maximum_power_mean_and_maybe_force_steinmetz(inputs, models);
    _coreLossesModelSteinmetz = CoreLossesModel::factory(models);
    _coreLossesModelProprietary = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Proprietary"}}));
    _magnetizingInductance = MagnetizingInductance(models["gapReluctance"]);
    _models = models;
}

std::pair<bool, double> MagneticFilterCoreDcAndSkinLosses::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    const auto& core = magnetic->get_core();

    if (inputs->get_operating_points().size() > 0 && magnetic->get_mutable_coil().get_functional_description().size() != inputs->get_operating_points()[0].get_excitations_per_winding().size()) {
        return {false, 0.0};
    }

    // This filter winds candidates with fast_wind() (no delimit/compact) for
    // speed. The flag lives on the process-global Settings singleton, so it must
    // be restored on every exit path — including the many throwing/early-return
    // paths below — or it leaks into all later code that reads it. RAII instead
    // of manual save/restore (which previously had no restore at all).
    SettingsGuard<bool> coilDelimitGuard(settings,
        &Settings::get_coil_delimit_and_compact,
        &Settings::set_coil_delimit_and_compact, false);

    std::string shapeName = core.get_shape_name();
    prepare_bobbin_for_non_pqi(magnetic, shapeName);

    auto currentNumberTurns = magnetic->get_coil().get_functional_description()[0].get_number_turns();
    NumberTurns numberTurns(currentNumberTurns);

    // Step size for the N sweep. With a ~10 iteration budget, a step of 1 only covers
    // N..N+10 (too narrow to find the loss minimum for larger designs). Using ~10% of
    // the starting N gives geometric-ish coverage out to ~2× N_start, spanning both
    // sides of the typical loss optimum.
    size_t numberTurnsStep = std::max<size_t>(1, static_cast<size_t>(std::ceil(currentNumberTurns * defaults.coreAdviserSkinEffectTurnsStepFactor)));

    std::vector<double> totalLossesPerOperatingPoint;
    std::vector<CoreLossesOutput> coreLossesPerOperatingPoint;
    std::vector<WindingLossesOutput> windingLossesPerOperatingPoint;
    double currentTotalLosses = std::numeric_limits<double>::infinity();
    double coreLosses = std::numeric_limits<double>::infinity();
    CoreLossesOutput coreLossesOutput;
    double ohmicAndSkinEffectLosses = std::numeric_limits<double>::infinity();
    WindingLossesOutput windingLossesOutput;
    windingLossesOutput.set_origin(ResultOrigin::SIMULATION);
    double newTotalLosses = std::numeric_limits<double>::infinity();
    auto previousNumberTurnsPrimary = currentNumberTurns;

    size_t iteration = defaults.coreAdviserSkinEffectMaxIterations;

    Coil coil = magnetic->get_coil();

    // Energy-storing inductors saturate as turns rise: at a fixed gap isat ∝ 1/N
    // (L ∝ N², isat = B_sat·N·A_e/L). The loss-minimising turn sweep below would
    // otherwise push N up (lower core loss) past the saturation-safe point — the
    // CoreAdviser's saturation-safe turn/gap sizing then gets silently undone and
    // the design saturates. Cap the sweep using the SAME
    // Magnetic::calculate_saturation_current the saturation filter and realism
    // gate use, so the loss-optimal turn count we return still clears the margin.
    bool isInductor = is_inductor(*inputs);
    double saturationMargin = Settings::GetInstance().get_core_adviser_saturation_margin();

    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
        auto operatingPoint = inputs->get_operating_point(operatingPointIndex);
        // Derating (ABT #13): hot junction corner, consistent with the saturation
        // filter's inductor gate (raw B_sat below).
        double temperature = saturation_derating_temperature(operatingPoint.get_conditions().get_ambient_temperature());
        OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];
        double saturationPeakCurrent = (isInductor && excitation.get_current() && excitation.get_current()->get_processed()
                                        && excitation.get_current()->get_processed()->get_peak())
                                       ? std::abs(excitation.get_current()->get_processed()->get_peak().value()) : 0.0;
        size_t numberTimeouts = 0;

        // Track the loss-minimum across the sweep. The do-while only gates progress;
        // the coil state at break time may be worse than an earlier iteration, so we
        // remember the best and restore it after the loop.
        double bestTotalLosses = std::numeric_limits<double>::infinity();
        uint64_t bestNumberTurnsPrimary = currentNumberTurns;
        CoreLossesOutput bestCoreLossesOutput;
        WindingLossesOutput bestWindingLossesOutput;
        bestWindingLossesOutput.set_origin(ResultOrigin::SIMULATION);

        do {
            currentTotalLosses = newTotalLosses;
            auto numberTurnsCombination = numberTurns.get_next_number_turns_combination(numberTurnsStep);
            coil.get_mutable_functional_description()[0].set_number_turns(numberTurnsCombination[0]);
            coil.fast_wind();  // delimit/compact disabled by coilDelimitGuard above

            // Saturation cap for energy-storing inductors: this (higher) turn count
            // is checked against the gap-aware saturation current; once it dips
            // below margin·I_pk, every higher N saturates too, so revert to the last
            // safe turn count and stop sweeping.
            if (saturationPeakCurrent > 0) {
                Magnetic saturationMagnetic = *magnetic;
                saturationMagnetic.set_coil(coil);
                double saturationCurrent = 0;
                try { saturationCurrent = saturationMagnetic.calculate_saturation_current(temperature, /*proportion=*/false); }
                catch (const std::exception& e) {
                    // ABT #121.4: was silent — with saturationCurrent = 0 the
                    // saturation cap on the turns sweep is disabled, so say so.
                    logEntry(std::string("Losses filter: saturation current failed during turns sweep: ") +
                             e.what() + " — saturation cap disabled for this sweep step", "CoreAdviser", 2);
                    saturationCurrent = 0;
                }
                if (saturationCurrent > 0 && saturationCurrent < saturationMargin * saturationPeakCurrent) {
                    coil.get_mutable_functional_description()[0].set_number_turns(previousNumberTurnsPrimary);
                    settings.set_coil_delimit_and_compact(false);
                    coil.fast_wind();
                    break;
                }
            }

            auto [magnetizingInductance, magneticFluxDensity] = _magnetizingInductance.calculate_inductance_and_magnetic_flux_density(core, coil, &operatingPoint);

            if (!check_requirement(inputs->get_design_requirements().get_magnetizing_inductance(), magnetizingInductance.get_magnetizing_inductance().get_nominal().value())) {
                if (resolve_dimensional_values(inputs->get_design_requirements().get_magnetizing_inductance()) < resolve_dimensional_values(magnetizingInductance.get_magnetizing_inductance().get_nominal().value())) {
                    coil.get_mutable_functional_description()[0].set_number_turns(previousNumberTurnsPrimary);
                    coil.fast_wind();  // delimit/compact disabled by coilDelimitGuard above
                    break;
                }
            }
            else {
                previousNumberTurnsPrimary = numberTurnsCombination[0];
            }

            if (!is_pqi_or_ui_shape(shapeName)) {
                if (!coil.get_turns_description()) {
                    // Phase 1 fix: previously silently broke with
                    // newTotalLosses = coreLosses (DBL_MAX initial),
                    // relying on later checks to drop the OP. fast_wind()
                    // failed to lay out the coil — reject explicitly.
                    return {false, 0.0};
                }
            }

            excitation.set_magnetic_flux_density(magneticFluxDensity);
            {
                auto pick = compute_core_losses_with_negative_guard(core, excitation, temperature, _coreLossesModelSteinmetz, _coreLossesModelProprietary);
                if (!pick.ok) {
                    // Phase 1 fix: was a silent `break`. Proprietary models
                    // can return tiny negatives (~µW) from interpolation
                    // noise at low excitation.
                    return {false, 0.0};
                }
                coreLossesOutput = pick.output;
                coreLosses = pick.value;
            }

            if (coreLosses < 0) {
                throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happend in core losses calculation for magnetic: " + magnetic->get_manufacturer_info().value().get_reference().value());
            }

            if (!coil.get_turns_description()) {
                // Phase 1 fix: explicit reject rather than silent break.
                return {false, 0.0};
            }

            if (!is_pqi_or_ui_shape(shapeName)) {
                windingLossesOutput = _windingOhmicLosses.calculate_ohmic_losses(coil, operatingPoint, temperature);
                windingLossesOutput = _windingSkinEffectLosses.calculate_skin_effect_losses(coil, temperature, windingLossesOutput, settings.get_harmonic_amplitude_threshold());

                ohmicAndSkinEffectLosses = windingLossesOutput.get_winding_losses();
                newTotalLosses = coreLosses + ohmicAndSkinEffectLosses;
                if (ohmicAndSkinEffectLosses < 0) {
                    throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happend in ohmic losses calculation for magnetic: " + magnetic->get_manufacturer_info().value().get_reference().value() + " ohmicAndSkinEffectLosses: " + std::to_string(ohmicAndSkinEffectLosses));
                }
            }
            else {
                // PQI / UI shapes: winding losses can't be computed without a
                // resolved coil turns_description here (these shapes use
                // integrated windings whose layout isn't produced by
                // fast_wind()). Use core-only losses for filtering and stop
                // sweeping turn counts — there's nothing to optimise against.
                // Not a silent fallback: explicit per-shape policy.
                newTotalLosses = coreLosses;
                break;
            }

            if (!std::isfinite(newTotalLosses)) {
                throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Too large losses");
            }

            // Record the best point seen so far in this operating-point sweep.
            if (newTotalLosses < bestTotalLosses) {
                bestTotalLosses = newTotalLosses;
                bestNumberTurnsPrimary = numberTurnsCombination[0];
                bestCoreLossesOutput = coreLossesOutput;
                bestWindingLossesOutput = windingLossesOutput;
            }

            iteration--;
            if (iteration <=0) {
                numberTimeouts++;
                break;
            }
        }
        while(newTotalLosses < currentTotalLosses * defaults.coreAdviserThresholdValidity);

        // Restore the best N from the sweep so downstream code sees the loss-optimal coil.
        if (std::isfinite(bestTotalLosses) &&
            coil.get_functional_description()[0].get_number_turns() != static_cast<double>(bestNumberTurnsPrimary)) {
            coil.get_mutable_functional_description()[0].set_number_turns(bestNumberTurnsPrimary);
            coil.fast_wind();  // delimit/compact disabled by coilDelimitGuard above
        }

        if (std::isfinite(bestTotalLosses)) {
            magnetic->set_coil(coil);
            totalLossesPerOperatingPoint.push_back(bestTotalLosses);
            coreLossesPerOperatingPoint.push_back(bestCoreLossesOutput);
            windingLossesPerOperatingPoint.push_back(bestWindingLossesOutput);
        }
        else if (std::isfinite(coreLosses) && coreLosses > 0) {
            // Fallback: no point ever became finite (e.g. PQI/UI shortcut path with only core losses).
            magnetic->set_coil(coil);
            currentTotalLosses = newTotalLosses;
            totalLossesPerOperatingPoint.push_back(currentTotalLosses);
            coreLossesPerOperatingPoint.push_back(coreLossesOutput);
            windingLossesPerOperatingPoint.push_back(windingLossesOutput);
        }
    }

    return finalize_losses_scoring(totalLossesPerOperatingPoint, coreLossesPerOperatingPoint, windingLossesPerOperatingPoint,
                                    magnetic, inputs, outputs, _maximumPowerMean);
}

MagneticFilterLosses::MagneticFilterLosses(std::map<std::string, std::string> models) {
    _models = models;
}

std::pair<bool, double> MagneticFilterLosses::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    double scoring = 0;

    if (inputs->get_operating_points().size() > 0 && magnetic->get_mutable_coil().get_functional_description().size() != inputs->get_operating_points()[0].get_excitations_per_winding().size()) {
        return {false, 0.0};
    }

    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
        auto operatingPoint = inputs->get_operating_points()[operatingPointIndex];
        auto temperature = operatingPoint.get_conditions().get_ambient_temperature();
        auto windingLosses = _magneticSimulator.calculate_winding_losses(operatingPoint, *magnetic, temperature);
        auto windingLossesValue = windingLosses.get_winding_losses();
        auto coreLosses = _magneticSimulator.calculate_core_losses(operatingPoint, *magnetic);
        auto coreLossesValue = coreLosses.get_core_losses();
        scoring += windingLossesValue + coreLossesValue;

        if (outputs != nullptr) {
            while (outputs->size() < operatingPointIndex + 1) {
                outputs->push_back(Outputs());
            }
            (*outputs)[operatingPointIndex].set_core_losses(coreLosses);
            (*outputs)[operatingPointIndex].set_winding_losses(windingLosses);
        }
    }

    scoring /= inputs->get_operating_points().size();

    return {true, scoring};
}

MagneticFilterLossesNoProximity::MagneticFilterLossesNoProximity(std::map<std::string, std::string> models) {
    _models = models;
}

std::pair<bool, double> MagneticFilterLossesNoProximity::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    double scoring = 0;

    if (inputs->get_operating_points().size() > 0 && magnetic->get_mutable_coil().get_functional_description().size() != inputs->get_operating_points()[0].get_excitations_per_winding().size()) {
        return {false, 0.0};
    }

    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
        auto operatingPoint = inputs->get_operating_points()[operatingPointIndex];
        auto temperature = operatingPoint.get_conditions().get_ambient_temperature();
        auto windingLosses = _windingOhmicLosses.calculate_ohmic_losses(magnetic->get_coil(), operatingPoint, temperature);
        // Phase 1 fix: was previously called twice on the same windingLosses
        // (double-counting the skin-effect contribution). One call is correct.
        windingLosses = _windingSkinEffectLosses.calculate_skin_effect_losses(magnetic->get_coil(), temperature, windingLosses, 0.5);
        // windingLosses = WindingLosses::combine_turn_losses(windingLosses, magnetic->get_coil());
        double windingLossesValue = windingLosses.get_winding_losses();

        auto coreLosses = _magneticSimulator.calculate_core_losses(operatingPoint, *magnetic);
        auto coreLossesValue = coreLosses.get_core_losses();
        scoring += windingLossesValue + coreLossesValue;

        if (outputs != nullptr) {
            while (outputs->size() < operatingPointIndex + 1) {
                outputs->push_back(Outputs());
            }
            (*outputs)[operatingPointIndex].set_core_losses(coreLosses);
            (*outputs)[operatingPointIndex].set_winding_losses(windingLosses);
        }
    }

    scoring /= inputs->get_operating_points().size();

    return {true, scoring};
}

} // namespace OpenMagnetics
