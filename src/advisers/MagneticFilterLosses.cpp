// Phase 5: loss-filter implementations extracted from MagneticFilter.cpp.
// All four classes (CoreAndDcLosses, CoreDcAndSkinLosses, Losses,
// LossesNoProximity) share the helpers in MagneticFilterInternal.h and
// nothing else, so they sit naturally in their own translation unit.
//
// Class declarations remain in advisers/MagneticFilter.h; this file
// supplies the definitions only.

#include "advisers/MagneticFilter.h"
#include "advisers/MagneticFilterInternal.h"
#include "constructive_models/Bobbin.h"
#include "constructive_models/NumberTurns.h"
#include "physical_models/ThermalNode.h"
#include "physical_models/WindingLosses.h"
#include "physical_models/WindingProximityEffectLosses.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "support/Exceptions.h"
#include "support/Utils.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <numbers>
#include <sstream>

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

std::pair<bool, double> MagneticFilterCoreDcAndSkinLosses::evaluate_stand_in_with_every_winding(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    SettingsGuard<bool> coilDelimitGuard(settings,
        &Settings::get_coil_delimit_and_compact,
        &Settings::set_coil_delimit_and_compact, false);

    // A coil lacking windings is the core adviser's one-winding stand-in. A strand that carries
    // its winding's current with copper to spare (18.6 mm at 50 Hz) is scored on the minimum
    // copper the window copper capacity screen counts (with_copper_sized_to_current without the
    // window fill: this scoring filter runs on every candidate, and the minimum copper scores
    // every candidate alike; the core-stage temperature gate fills the window).
    Magnetic completed = with_copper_sized_to_current(with_every_winding(*magnetic, *inputs), *inputs,
                                                      defaults.maximumEffectiveCurrentDensity, false);
    const auto& core = completed.get_core();
    const std::string shapeName = core.get_shape_name();
    prepare_bobbin_for_non_pqi(&completed, shapeName);
    auto& coil = completed.get_mutable_coil();
    const bool windable = !is_pqi_or_ui_shape(shapeName);
    if (windable) {
        coil.fast_wind();  // delimit/compact disabled by coilDelimitGuard above
        if (!coil.get_turns_description()) {
            // As for one winding: a coil fast_wind() cannot lay out is rejected explicitly.
            return {false, 0.0};
        }
    }

    std::vector<double> totalLossesPerOperatingPoint;
    std::vector<CoreLossesOutput> coreLossesPerOperatingPoint;
    std::vector<WindingLossesOutput> windingLossesPerOperatingPoint;
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
        auto operatingPoint = inputs->get_operating_point(operatingPointIndex);
        // The same hot corner the one-winding path evaluates at.
        double temperature = saturation_derating_temperature(operatingPoint.get_conditions().get_ambient_temperature());
        auto [magnetizingInductance, magneticFluxDensity] = _magnetizingInductance.calculate_inductance_and_magnetic_flux_density(core, coil, &operatingPoint);
        (void) magnetizingInductance;
        OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];
        excitation.set_magnetic_flux_density(magneticFluxDensity);
        auto pick = compute_core_losses_with_negative_guard(core, excitation, temperature, _coreLossesModelSteinmetz, _coreLossesModelProprietary);
        if (!pick.ok) {
            return {false, 0.0};
        }
        if (pick.value < 0) {
            throw CalculationException(ErrorCode::CALCULATION_ERROR, "Negative core losses for magnetic: " + magnetic->get_reference());
        }
        WindingLossesOutput windingLossesOutput;
        windingLossesOutput.set_origin(ResultOrigin::SIMULATION);
        double totalLosses = pick.value;
        if (windable) {
            // PQI / UI shapes: the same explicit per-shape policy as the one-winding path
            // (core losses only, their integrated windings are not laid out by fast_wind()).
            windingLossesOutput = _windingOhmicLosses.calculate_ohmic_losses(coil, operatingPoint, temperature);
            windingLossesOutput = _windingSkinEffectLosses.calculate_skin_effect_losses(coil, temperature, windingLossesOutput, settings.get_harmonic_amplitude_threshold());
            if (windingLossesOutput.get_winding_losses() < 0) {
                throw CalculationException(ErrorCode::CALCULATION_ERROR, "Negative winding losses for magnetic: " + magnetic->get_reference());
            }
            totalLosses += windingLossesOutput.get_winding_losses();
        }
        if (!std::isfinite(totalLosses)) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Too large losses");
        }
        totalLossesPerOperatingPoint.push_back(totalLosses);
        coreLossesPerOperatingPoint.push_back(pick.output);
        windingLossesPerOperatingPoint.push_back(windingLossesOutput);
    }
    return finalize_losses_scoring(totalLossesPerOperatingPoint, coreLossesPerOperatingPoint, windingLossesPerOperatingPoint,
                                   magnetic, inputs, outputs, _maximumPowerMean);
}

std::pair<bool, double> MagneticFilterCoreDcAndSkinLosses::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    const auto& core = magnetic->get_core();

    if (inputs->get_operating_points().size() > 0 && magnetic->get_mutable_coil().get_functional_description().size() != inputs->get_operating_points()[0].get_excitations_per_winding().size()) {
        if (magnetic->get_coil().get_functional_description().size() == 1) {
            return evaluate_stand_in_with_every_winding(magnetic, inputs, outputs);
        }
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

// ---------------------------------------------------------------------------
// ABT #1426: loss-optimal (turns, gap) of a gapped inductor
// ---------------------------------------------------------------------------
//
// The CoreAdviser seeds an inductor at the FEWEST turns that hold its inductance and clear
// the saturation margin (add_initial_turns_by_inductance, ABT #13). For a DC-biased choke
// that is usually where the losses are lowest too: its flux swing is small, core loss is
// small, and every extra turn only adds copper. For an inductor whose AC flux swing
// dominates (an LLC/CLLC resonant inductor, a boost inductor in DCM, ...) it is the worst
// place to stop: B sits at the saturation ceiling and the core dissipates 4-10x what the
// copper does (ABT #1426: 9 turns, 0.34 T, 60 W core against 8 W copper).
//
// At a FIXED inductance L the two losses move in opposite directions with the turn count:
//
//   B(t) = L i(t) / (N Ae)                  -> core loss  ~ N^-beta   (beta ~ 2-3, Steinmetz)
//   window-filled copper per turn = k Aw / N -> copper loss ~ N^2
//
// so the total has a single minimum, at Pcore = (2 / beta) Pcu for a pure power law. The
// gap is what keeps L fixed as N grows (reluctance N^2 / L), so the pair is chosen together:
// the turn count by the loss estimate, then the gap re-solved to hold L at that count.
//
// The search never goes below the seeded turn count (the saturation floor): more turns at
// the same inductance only LOWER the peak flux density, so every candidate it considers is
// at least as far from saturation as the seed. It stops where the design stops being
// buildable: the re-solved gap must fit the column, stay under the post-gap fringing
// ceiling (single gap, or a distributed gap when the settings allow one), keep the
// inductance in its band, pass the saturation gate, and the adviser's stand-in winding must
// still fit the window.
//
// Cost: the loss estimate is analytical (no winding, no simulation). Core loss at L fixed
// does not depend on the gap, so the whole turn search runs WITHOUT gap solves; the gap is
// solved only for the chosen count (and a few bisection steps when that count is not
// buildable).
MagneticFilterInductorTurnsAndGapByLosses::MagneticFilterInductorTurnsAndGapByLosses(Inputs inputs, std::map<std::string, std::string> models)
    : _fringingFactorFilter(inputs, models) {
    _models = models;
    _magnetizingInductance = MagnetizingInductance(models["gapReluctance"]);
    _coreLossesModelSteinmetz = CoreLossesModel::factory(models);
    _coreLossesModelProprietary = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Proprietary"}}));
    // The same post-gap fringing ceiling reject_winding_killing_gaps applies to every
    // finalised gap: a re-solved gap past it would be culled right after this step.
    _fringingFactorFilter.set_fringing_factor_limit(defaults.coreAdviserWindingKillingFringingFactorLimit);

    _targetInductance = resolve_dimensional_values(inputs.get_design_requirements().get_magnetizing_inductance());
    if (!(_targetInductance > 0)) {
        throw InvalidInputException(ErrorCode::MISSING_DATA,
            "Loss-optimal inductor turns: the design has no positive magnetizing inductance target");
    }

    // Turns ratios of the other windings (primary / winding w). The copper of every winding
    // shares the one window, so the window-filled copper estimate below needs the
    // primary-referred ampere-turns of all of them.
    std::vector<double> turnsRatios;
    for (const auto& turnsRatio : inputs.get_design_requirements().get_turns_ratios()) {
        turnsRatios.push_back(resolve_dimensional_values(turnsRatio));
    }
    _turnsRatios = turnsRatios;
    _windingCurrents.resize(turnsRatios.size() + 1);

    const auto& operatingPoints = inputs.get_operating_points();
    if (operatingPoints.empty()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA,
            "Loss-optimal inductor turns: the design has no operating points");
    }
    for (size_t operatingPointIndex = 0; operatingPointIndex < operatingPoints.size(); ++operatingPointIndex) {
        const auto& operatingPoint = operatingPoints[operatingPointIndex];
        // The same temperature MagneticFilterCoreDcAndSkinLosses evaluates both losses at
        // (the hot corner, ABT #13), so this choice and the loss ranking that follows it
        // judge the design at the same point.
        _temperatures.push_back(saturation_derating_temperature(operatingPoint.get_conditions().get_ambient_temperature()));

        // Optimal window allocation (each winding gets the share of copper proportional to
        // its ampere-turns) gives a total copper loss of
        //     P = rho MLT / (k Aw) (sum_w N_w I_w,rms)^2 = rho MLT N^2 / (k Aw) (sum_w I_w,rms / n_w)^2
        // with n_w = N / N_w. For a single winding this is just N^2 I_rms^2 rho MLT / (k Aw).
        const auto& excitations = operatingPoint.get_excitations_per_winding();
        if (excitations.size() > turnsRatios.size() + 1) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "Loss-optimal inductor turns: operating point " + std::to_string(operatingPointIndex) +
                " has " + std::to_string(excitations.size()) + " excitations but only " +
                std::to_string(turnsRatios.size()) + " turns ratios");
        }
        double primaryReferredCurrentRms = 0;
        for (size_t windingIndex = 0; windingIndex < excitations.size(); ++windingIndex) {
            const auto& excitation = excitations[windingIndex];
            if (!excitation.get_current() || !excitation.get_current()->get_processed() ||
                !excitation.get_current()->get_processed()->get_rms()) {
                throw InvalidInputException(ErrorCode::MISSING_DATA,
                    "Loss-optimal inductor turns: operating point " + std::to_string(operatingPointIndex) +
                    ", winding " + std::to_string(windingIndex) + " has no processed RMS current");
            }
            _windingCurrents[windingIndex].push_back(excitation.get_current().value());
            double currentRms = excitation.get_current()->get_processed()->get_rms().value();
            double turnsRatio = windingIndex == 0 ? 1.0 : turnsRatios[windingIndex - 1];
            primaryReferredCurrentRms += currentRms / turnsRatio;
        }
        _primaryReferredCurrentsRms.push_back(primaryReferredCurrentRms);

        // Harmonics of the window's ampere-turns for the proximity estimate: the primary
        // current's spectrum, scaled to the primary-referred total above (for one winding
        // the scale is 1). DC carries no proximity loss; harmonics below the settings'
        // amplitude threshold are dropped as every winding-loss path in MKF drops them.
        auto primaryCurrent = excitations[0].get_current().value();
        if (!primaryCurrent.get_harmonics() || !primaryCurrent.get_processed()->get_effective_frequency()) {
            throw InvalidInputException(ErrorCode::MISSING_DATA,
                "Loss-optimal inductor turns: operating point " + std::to_string(operatingPointIndex) +
                ", primary current has no harmonics or effective frequency");
        }
        _maximumEffectiveFrequency = std::max(_maximumEffectiveFrequency, primaryCurrent.get_processed()->get_effective_frequency().value());
        double primaryCurrentRms = primaryCurrent.get_processed()->get_rms().value();
        auto harmonics = primaryCurrent.get_harmonics().value();
        const auto& amplitudes = harmonics.get_amplitudes();
        const auto& frequencies = harmonics.get_frequencies();
        double maximumAmplitude = 0;
        for (size_t harmonicIndex = 1; harmonicIndex < amplitudes.size(); ++harmonicIndex) {
            maximumAmplitude = std::max(maximumAmplitude, amplitudes[harmonicIndex]);
        }
        std::vector<double> keptAmplitudes;
        std::vector<double> keptFrequencies;
        if (primaryCurrentRms > 0) {
            double scale = primaryReferredCurrentRms / primaryCurrentRms;
            for (size_t harmonicIndex = 1; harmonicIndex < amplitudes.size(); ++harmonicIndex) {
                if (amplitudes[harmonicIndex] < maximumAmplitude * settings.get_harmonic_amplitude_threshold()) {
                    continue;
                }
                keptAmplitudes.push_back(amplitudes[harmonicIndex] * scale);
                keptFrequencies.push_back(frequencies[harmonicIndex]);
            }
        }
        _primaryReferredCurrentHarmonicAmplitudes.push_back(keptAmplitudes);
        _primaryReferredCurrentHarmonicFrequencies.push_back(keptFrequencies);
    }
}

bool MagneticFilterInductorTurnsAndGapByLosses::applies_to_design(const Inputs& inputs) {
    if (!is_inductor(inputs) || inputs.get_operating_points().empty()) {
        return false;
    }
    // Suppression chokes are sized by impedance, not by an inductance target and losses.
    if (is_tagged_common_mode_choke(inputs, inputs.get_operating_points()[0])) {
        return false;
    }
    auto application = inputs.get_design_requirements().get_application();
    if (application && application.value() == "interferenceSuppression") {
        return false;
    }
    const auto& inductance = inputs.get_design_requirements().get_magnetizing_inductance();
    bool hasTarget = (inductance.get_nominal() && inductance.get_nominal().value() > 0) ||
                     (inductance.get_maximum() && inductance.get_maximum().value() > 0) ||
                     (inductance.get_minimum() && inductance.get_minimum().value() > 0);
    return hasTarget;
}

bool MagneticFilterInductorTurnsAndGapByLosses::applies_to_candidate(const Magnetic& magnetic) {
    const auto& core = magnetic.get_core();
    // A toroid (powder) has no discrete gap to re-solve; PQI / UI integrated windings are
    // not laid out by fast_wind, so the stand-in winding's fit cannot be checked on them.
    if (core.get_shape_family() == CoreShapeFamily::T || is_pqi_or_ui_shape(core.get_shape_name())) {
        return false;
    }
    if (magnetic.get_coil().get_functional_description().size() != 1) {
        return false;
    }
    return true;
}

double MagneticFilterInductorTurnsAndGapByLosses::calculate_core_losses(const Core& core, double numberTurns,
                                                                        std::vector<OperatingPoint>& preparedOperatingPoints) {
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    double coreLosses = 0;
    for (size_t operatingPointIndex = 0; operatingPointIndex < preparedOperatingPoints.size(); ++operatingPointIndex) {
        auto& operatingPoint = preparedOperatingPoints[operatingPointIndex];
        // At the target inductance the core carries flux L i / N whatever the gap, which is
        // exactly what the gap is re-solved to deliver at this turn count.
        auto magneticFluxDensity = MagnetizingInductance::calculate_magnetic_flux_density_from_inductance(
            _targetInductance, numberTurns, effectiveArea, &operatingPoint);
        OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];
        excitation.set_magnetic_flux_density(magneticFluxDensity);
        auto pick = compute_core_losses_with_negative_guard(core, excitation, _temperatures[operatingPointIndex],
                                                           _coreLossesModelSteinmetz, _coreLossesModelProprietary);
        if (!pick.ok) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        if (!std::isfinite(pick.value)) {
            throw CalculationException(ErrorCode::CALCULATION_ERROR,
                "Loss-optimal inductor turns: core loss estimate of " + core.get_name().value_or("?") + " at " +
                std::to_string(static_cast<int64_t>(numberTurns)) + " turns, operating point " +
                std::to_string(operatingPointIndex) + " is not finite (" + std::to_string(pick.value) + " W)");
        }
        if (pick.value < 0) {
            throw CalculationException(ErrorCode::CALCULATION_ERROR,
                "Loss-optimal inductor turns: negative core losses for core " + core.get_name().value_or("?"));
        }
        coreLosses += pick.value;
    }
    // Operating points weigh equally, as in finalize_losses_scoring.
    return coreLosses / preparedOperatingPoints.size();
}

std::pair<bool, double> MagneticFilterInductorTurnsAndGapByLosses::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, [[maybe_unused]] std::vector<Outputs>* outputs) {
    if (!applies_to_candidate(*magnetic)) {
        return {true, 0};
    }
    Core core = magnetic->get_core();
    if (!core.get_processed_description()) {
        throw CoreNotProcessedException("Loss-optimal inductor turns: core " + core.get_name().value_or("?") + " is not processed");
    }
    Coil coil = magnetic->get_coil();
    auto seededNumberTurns = static_cast<int64_t>(std::llround(coil.get_functional_description()[0].get_number_turns()));
    if (seededNumberTurns < 1) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "Loss-optimal inductor turns: candidate " + core.get_name().value_or("?") + " has no seeded turns");
    }

    // The operating points with the magnetizing current the flux follows. For one winding that
    // winding's current IS the magnetizing current, and computing the seed's inductance and flux
    // once puts it there in the form the inductance model uses (standardised, power-of-two).
    // For several windings it is not: the stand-in coil has one winding, and the inductance
    // model would take the primary's own current as the magnetizing current. A flyback primary
    // carries only its on-time share and drops to zero when the secondary takes over, a step
    // the flux never makes; fed to the core-loss models as B it was a flux discontinuity, and
    // the iGSE integral returned NaN on it. The design's magnetizing current, derived from the
    // primary volt-seconds at the target inductance (pre_process_inputs), is the right one.
    std::vector<OperatingPoint> preparedOperatingPoints;
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
        auto operatingPoint = inputs->get_operating_point(operatingPointIndex);
        if (operatingPoint.get_excitations_per_winding().size() > 1) {
            if (!operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current()) {
                throw InvalidInputException(ErrorCode::MISSING_DATA,
                    "Loss-optimal inductor turns: operating point " + std::to_string(operatingPointIndex) +
                    " has several windings but no magnetizing current (the adviser derives it in pre_process_inputs)");
            }
        }
        else {
            _magnetizingInductance.calculate_inductance_and_magnetic_flux_density(core, coil, &operatingPoint);
        }
        preparedOperatingPoints.push_back(operatingPoint);
    }

    // Copper: the window filled with litz. At the core stage no wire is chosen yet, so the
    // estimate is the copper a winding adviser can put in this window at best:
    //   - strands of diameter delta / 3 at the highest effective frequency -- the middle of
    //     the delta/2 .. delta/4 grid MKF's WireAdviser synthesises litz from
    //     (synthesize_litz_candidates, Sullivan); skin effect inside such a strand is small;
    //   - packed at k = Wire::get_filling_factor_round(strand) x Bobbin::get_filling_factor
    //     (the window utilisation MagneticFilterAreaProduct uses), so each of N turns gets
    //     k Aw / N of copper: DC loss = rho MLT N^2 I_rms^2 / (k Aw);
    //   - MLT = the turn length at the middle of the bobbin window (the turn-length geometry
    //     the coil winder charges, ThermalNetworkNode::computeTurnLengthAtRadius);
    //   - proximity from the one-dimensional window field of a window-filled winding: the
    //     field grows linearly across the winding build to N I / b at its outer face (b the
    //     winding breadth along the column, Dowell's assumption), so its mean square over
    //     the strands is (N I_h)^2 / (3 b^2); each strand loses what MKF's litz proximity
    //     model gives at that field. Strand count x per-strand loss is also ~N^2, because the
    //     strand count at a fixed fill does not depend on N.
    // Both terms scale as N^2, so the copper is one coefficient per candidate times N^2.
    // Gap fringing near the winding is not in it (the post-gap fringing ceiling bounds it).
    std::string shapeName = core.get_shape_name();
    Magnetic magneticWithBobbin = *magnetic;
    prepare_bobbin_for_non_pqi(&magneticWithBobbin, shapeName);
    auto bobbin = magneticWithBobbin.get_mutable_coil().resolve_bobbin();
    auto bobbinProcessed = bobbin.get_processed_description().value();
    auto bobbinWindow = bobbinProcessed.get_winding_windows()[0];
    if (!bobbinProcessed.get_column_width() || !bobbinWindow.get_coordinates() || !bobbinWindow.get_height()) {
        throw CoilNotProcessedException("Loss-optimal inductor turns: quick bobbin of " + core.get_name().value_or("?") +
                                        " has no column width, window coordinates or window height");
    }
    double meanTurnLength = ThermalNetworkNode::computeTurnLengthAtRadius(
        bobbinWindow.get_coordinates().value()[0], bobbinProcessed.get_column_shape(),
        bobbinProcessed.get_column_width().value(), bobbinProcessed.get_column_depth());
    double windingBreadth = bobbinWindow.get_height().value();
    auto coreWindow = core.get_winding_windows()[0];
    double coreWindowArea = coreWindow.get_area().value();
    double bobbinFillingFactor = Bobbin::get_filling_factor(coreWindow.get_width().value(), coreWindow.get_height().value());
    auto wire = coil.get_wires()[0];

    double hottestTemperature = *std::max_element(_temperatures.begin(), _temperatures.end());
    double strandDiameter = WindingSkinEffectLosses::calculate_skin_depth(wire, _maximumEffectiveFrequency, hottestTemperature) / 3;
    auto strand = Wire::create_quick_litz_wire(strandDiameter, 1);
    double litzCopperArea = Wire::get_filling_factor_round(strandDiameter) * bobbinFillingFactor * coreWindowArea;
    double strandConductingArea = std::numbers::pi * pow(strandDiameter, 2) / 4;
    double numberStrands = litzCopperArea / strandConductingArea;
    auto proximityModel = WindingProximityEffectLosses::get_model(WireType::LITZ);

    double meanCopperLossesPerTurnSquared = 0;
    for (size_t operatingPointIndex = 0; operatingPointIndex < _temperatures.size(); ++operatingPointIndex) {
        double temperature = _temperatures[operatingPointIndex];
        double resistivity = WindingOhmicLosses::calculate_dc_resistance_per_meter(wire, temperature) *
                             wire.calculate_conducting_area();
        double dcLosses = pow(_primaryReferredCurrentsRms[operatingPointIndex], 2) * resistivity * meanTurnLength / litzCopperArea;
        double proximityLosses = 0;
        const auto& amplitudes = _primaryReferredCurrentHarmonicAmplitudes[operatingPointIndex];
        const auto& frequencies = _primaryReferredCurrentHarmonicFrequencies[operatingPointIndex];
        for (size_t harmonicIndex = 0; harmonicIndex < amplitudes.size(); ++harmonicIndex) {
            ComplexFieldPoint fieldPoint;
            // Root-mean-square over the winding build of the field at ONE turn (N = 1).
            fieldPoint.set_real(amplitudes[harmonicIndex] / (windingBreadth * sqrt(3)));
            fieldPoint.set_imaginary(0);
            fieldPoint.set_point(std::vector<double>{0, 0});
            double strandLossesPerMeter = proximityModel->calculate_turn_losses(strand, frequencies[harmonicIndex], {fieldPoint}, temperature);
            proximityLosses += strandLossesPerMeter * meanTurnLength * numberStrands;
        }
        meanCopperLossesPerTurnSquared += dcLosses + proximityLosses;
    }
    meanCopperLossesPerTurnSquared /= _temperatures.size();
    if (!std::isfinite(meanCopperLossesPerTurnSquared) || meanCopperLossesPerTurnSquared < 0) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR,
            "Loss-optimal inductor turns: copper loss estimate of " + core.get_name().value_or("?") +
            " is not a finite non-negative number (" + std::to_string(meanCopperLossesPerTurnSquared) + " W per turn squared)");
    }

    // Most turns the window can carry. The adviser's stand-in winding (get_dummy_coil) is what
    // the loss ranking winds next and what the coil adviser has to beat: strands two skin
    // depths across, as many in parallel as keep every winding under both MKF current-density
    // ceilings (maximumEffectiveCurrentDensity and maximumCurrentDensity, the same sizing
    // get_dummy_coil applies), at least one. Every winding shares the window, winding w
    // carrying N / n_w turns, so at N primary turns the conductors need
    //     N sum_w (parallels_w A_strand / n_w)
    // of copper, packed at the stand-in strand's filling factor inside the usable bobbin
    // fraction. That bounds N for multi-winding designs too, where the one-winding stand-in
    // coil cannot be wound to check it.
    double standInConductingDiameter = resolve_dimensional_values(wire.get_conducting_diameter().value());
    double standInCopperArea = Wire::get_filling_factor_round(standInConductingDiameter) * bobbinFillingFactor * coreWindowArea;
    double standInStrandArea = wire.calculate_conducting_area();
    double conductorAreaPerPrimaryTurn = 0;
    for (size_t windingIndex = 0; windingIndex < _windingCurrents.size(); ++windingIndex) {
        int parallels = 1;
        for (auto& current : _windingCurrents[windingIndex]) {
            int parallelsEffective = Wire::calculate_number_parallels_needed(current, hottestTemperature, wire, defaults.maximumEffectiveCurrentDensity);
            double dcCurrentDensity = wire.calculate_dc_current_density(current);
            int parallelsDc = static_cast<int>(std::ceil(dcCurrentDensity / defaults.maximumCurrentDensity));
            parallels = std::max({parallels, parallelsEffective, parallelsDc});
        }
        double turnsRatio = windingIndex == 0 ? 1.0 : _turnsRatios[windingIndex - 1];
        conductorAreaPerPrimaryTurn += parallels * standInStrandArea / turnsRatio;
    }
    auto maximumNumberTurnsByWindow = static_cast<int64_t>(std::floor(standInCopperArea / conductorAreaPerPrimaryTurn));
    if (maximumNumberTurnsByWindow <= seededNumberTurns) {
        return {true, 0};
    }

    bool hasTurnsRatios = !inputs->get_design_requirements().get_turns_ratios().empty();
    auto snapToTurnsRatios = [&](int64_t numberTurns) -> int64_t {
        if (!hasTurnsRatios) {
            return numberTurns;
        }
        NumberTurns numberTurnsCombination(numberTurns, inputs->get_design_requirements());
        return static_cast<int64_t>(numberTurnsCombination.get_next_number_turns_combination()[0]);
    };

    std::map<int64_t, double> totalLossesCache;
    bool unscorable = false;
    auto totalLosses = [&](int64_t numberTurns) -> double {
        numberTurns = snapToTurnsRatios(numberTurns);
        auto it = totalLossesCache.find(numberTurns);
        if (it != totalLossesCache.end()) {
            return it->second;
        }
        double coreLosses = calculate_core_losses(core, numberTurns, preparedOperatingPoints);
        if (std::isnan(coreLosses)) {
            unscorable = true;
            coreLosses = std::numeric_limits<double>::infinity();
        }
        double losses = coreLosses + meanCopperLossesPerTurnSquared * pow(numberTurns, 2);
        totalLossesCache[numberTurns] = losses;
        return losses;
    };

    // Bracket the minimum by geometric steps up from the seed (each step +25 %: a 3x range
    // is bracketed in ~5 evaluations), stopping as soon as the total loss rises.
    const double bracketGrowth = 1.25;
    int64_t lower = seededNumberTurns;
    int64_t best = seededNumberTurns;
    double bestLosses = totalLosses(seededNumberTurns);
    if (unscorable) {
        logEntry("Loss-optimal inductor turns: core losses of " + core.get_name().value_or("?") +
                 " cannot be computed (no Steinmetz or proprietary data); keeping the seeded turns, "
                 "the loss ranking will judge it", "CoreAdviser", 2);
        return {true, 0};
    }
    int64_t upper = best;
    while (true) {
        int64_t next = std::min(maximumNumberTurnsByWindow,
                                std::max(best + 1, static_cast<int64_t>(std::ceil(best * bracketGrowth))));
        if (next <= best) {
            upper = best;
            break;
        }
        double nextLosses = totalLosses(next);
        if (nextLosses >= bestLosses) {
            upper = next;
            break;
        }
        lower = best;
        best = next;
        bestLosses = nextLosses;
        if (next == maximumNumberTurnsByWindow) {
            upper = next;
            break;
        }
    }
    // Refine inside [lower, upper] by integer ternary search (the total is unimodal in N).
    while (upper - lower > 2) {
        int64_t third = (upper - lower) / 3;   // integer index arithmetic, not physics
        int64_t left = lower + third;
        int64_t right = upper - third;
        if (totalLosses(left) < totalLosses(right)) {
            upper = right;
        }
        else {
            lower = left;
        }
    }
    for (int64_t numberTurns = lower; numberTurns <= upper; ++numberTurns) {
        double losses = totalLosses(numberTurns);
        if (losses < bestLosses) {
            best = numberTurns;
            bestLosses = losses;
        }
    }
    best = snapToTurnsRatios(best);
    if (best <= seededNumberTurns) {
        return {true, bestLosses};
    }

    // Build the design at a turn count: gap re-solved to hold the inductance, single gap
    // first, distributed when the single one does not fit or fringes past the ceiling and
    // the settings allow it. nullopt when no gap makes it buildable.
    Inputs gappingInputs;
    gappingInputs.set_design_requirements(inputs->get_design_requirements());
    std::vector<GappingType> gappingTypes{GappingType::GROUND};
    if (settings.get_core_adviser_include_distributed_gaps()) {
        gappingTypes.push_back(GappingType::DISTRIBUTED);
    }
    // A single gap that cannot be built (does not fit the column, or fringes past the
    // ceiling) at N turns cannot be built at more turns either -- the gap only grows -- so
    // once one fails the single gap is not solved again above that count.
    int64_t singleGapUnbuildableFrom = std::numeric_limits<int64_t>::max();
    auto buildDesign = [&](int64_t numberTurns) -> std::optional<std::pair<Magnetic, GappingType>> {
        Coil trialCoil = coil;
        trialCoil.get_mutable_functional_description()[0].set_number_turns(numberTurns);
        for (auto gappingType : gappingTypes) {
            bool singleGap = gappingType == GappingType::GROUND;
            if (singleGap && numberTurns >= singleGapUnbuildableFrom) {
                continue;
            }
            auto markUnbuildable = [&]() {
                if (singleGap) {
                    singleGapUnbuildableFrom = std::min(singleGapUnbuildableFrom, numberTurns);
                }
            };
            std::vector<CoreGap> gaps;
            try {
                gaps = _magnetizingInductance.calculate_gapping_from_number_turns_and_inductance(core, trialCoil, &gappingInputs, gappingType);
            }
            catch (const GapException&) {
                // The gap this turn count needs cannot be built on this core: a verdict on
                // this (N, gap), not an engine error.
                markUnbuildable();
                continue;
            }
            if (gaps.empty()) {
                markUnbuildable();
                continue;
            }
            if (gappingType == GappingType::DISTRIBUTED) {
                // A distributed gap solved down to the residual length is no gap to machine:
                // the extra reluctance would come only from counting contact faces.
                bool residualOnly = true;
                for (const auto& gap : gaps) {
                    if (gap.get_type() != GapType::RESIDUAL && gap.get_length() > Constants().residualGap) {
                        residualOnly = false;
                    }
                }
                if (residualOnly) {
                    continue;
                }
            }
            Core trialCore = core;
            trialCore.set_gapping(gaps);
            if (!trialCore.process_gap()) {
                markUnbuildable();
                continue;
            }
            Magnetic trial = *magnetic;
            trial.set_core(trialCore);
            trial.set_coil(trialCoil);
            if (!_fringingFactorFilter.evaluate_magnetic(&trial, inputs).first) {
                markUnbuildable();
                continue;
            }
            if (!_inductanceFilter.evaluate_magnetic(&trial, inputs).first) {
                continue;
            }
            if (!_saturationFilter.evaluate_magnetic(&trial, inputs).first) {
                continue;
            }
            // The stand-in winding must still fit: the loss ranking winds it next and drops
            // the candidate if it does not.
            {
                SettingsGuard<bool> coilDelimitGuard(settings,
                    &Settings::get_coil_delimit_and_compact,
                    &Settings::set_coil_delimit_and_compact, false);
                Magnetic wound = trial;
                prepare_bobbin_for_non_pqi(&wound, shapeName);
                auto woundCoil = wound.get_coil();
                woundCoil.fast_wind();
                if (!woundCoil.get_turns_description()) {
                    // The window is full at this count, whatever the gap.
                    return std::nullopt;
                }
            }
            return std::pair<Magnetic, GappingType>{trial, gappingType};
        }
        return std::nullopt;
    };

    // Buildability only gets harder with more turns (a longer gap), so find the largest
    // buildable count in (seed, best] by bisection. Losses fall monotonically from the seed
    // up to `best`, so that count is the constrained optimum.
    std::optional<std::pair<Magnetic, GappingType>> chosen = buildDesign(best);
    int64_t chosenNumberTurns = best;
    if (!chosen) {
        int64_t buildable = seededNumberTurns;
        int64_t unbuildable = best;
        while (unbuildable - buildable > 1) {
            int64_t middle = snapToTurnsRatios(buildable + (unbuildable - buildable) / 2);   // index arithmetic
            if (middle >= unbuildable || middle <= buildable) {
                break;
            }
            auto trial = buildDesign(middle);
            if (trial) {
                buildable = middle;
                chosen = trial;
                chosenNumberTurns = middle;
            }
            else {
                unbuildable = middle;
            }
        }
    }
    if (!chosen) {
        return {true, bestLosses};
    }

    // Name the new gap the way add_gapping_standard_cores names it.
    Core chosenCore = chosen->first.get_core();
    auto name = chosenCore.get_name().value_or("unnamed");
    size_t gappedPosition = name.find(" gapped ");
    if (gappedPosition != std::string::npos) {
        name = name.substr(0, gappedPosition);
    }
    size_t ungappedPosition = name.find(" ungapped");
    if (ungappedPosition != std::string::npos) {
        name = name.substr(0, ungappedPosition);
    }
    // The gaps the solver cut (ground or distributed: SUBTRACTIVE), not the residual
    // contact gaps of the other columns.
    std::vector<CoreGap> cutGaps;
    for (const auto& gap : chosenCore.get_gapping()) {
        if (gap.get_type() != GapType::RESIDUAL) {
            cutGaps.push_back(gap);
        }
    }
    if (cutGaps.empty()) {
        throw GapException("Loss-optimal inductor turns: the re-solved gapping of " + name + " has no cut gap");
    }
    std::stringstream gapDescription;
    gapDescription << std::fixed << std::setprecision(2);
    if (chosen->second == GappingType::DISTRIBUTED) {
        gapDescription << cutGaps.size() << " x ";
    }
    gapDescription << cutGaps[0].get_length() * 1000;
    chosenCore.set_name(name + " gapped " + gapDescription.str() + " mm");
    magnetic->set_core(chosenCore);
    magnetic->get_mutable_coil().get_mutable_functional_description()[0].set_number_turns(chosenNumberTurns);
    logEntry("Loss-optimal inductor turns: " + chosenCore.get_name().value_or("?") + " from " +
             std::to_string(seededNumberTurns) + " to " + std::to_string(chosenNumberTurns) + " turns",
             "CoreAdviser", 2);
    return {true, totalLosses(chosenNumberTurns)};
}

} // namespace OpenMagnetics
