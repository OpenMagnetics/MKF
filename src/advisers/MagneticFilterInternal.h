// Phase 5: private helpers shared between MagneticFilter.cpp and any
// future split-out filter implementation files. Header-only / inline so
// each .cpp gets its own copy without ODR collisions, and so static-
// initialisation order vs. the `defaults` / `settings` globals in
// support/Utils.h is well-defined (those are inline globals, fine to
// reference from inline helpers).
//
// Nothing here is part of the public MagneticFilter API; do not include
// from outside src/advisers/.
#pragma once

#include "advisers/MagneticFilter.h"
#include "constructive_models/Bobbin.h"
#include "physical_models/CoreLosses.h"
#include "physical_models/MagnetizingInductance.h"
#include "processors/Inputs.h"
#include "Defaults.h"
#include "support/Exceptions.h"
#include "support/Utils.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <optional>
#include <utility>
#include <vector>

namespace OpenMagnetics {

// The role of the magnetic (Inputs.h, magnetic_role) decides only the categorical choices: gapping,
// the inductor loss optimiser, the ferrite-toroid drop and the choke handling. The magnitudes
// (saturation, the turn seed, the energy check, the loss-sweep cap) read the magnetizing current:
// on a single-winding design the winding current, on any other the one the magnetizing-inductance
// model stamps from the circuit (Faraday on the voltage, or the ampere-turns of the currents).
inline bool is_transformer_role(const Inputs& inputs) {
    return magnetic_role(inputs).role == MagneticRole::TRANSFORMER;
}

// The peak of the magnetizing current stamped on an operating point's first excitation (by the
// magnetizing-inductance model or Inputs::process): the largest magnitude of its waveform, else the
// magnitude of its processed peak. Throws when there is none.
inline double stamped_magnetizing_current_peak(const OperatingPoint& operatingPoint) {
    const auto& excitations = operatingPoint.get_excitations_per_winding();
    if (excitations.empty() || !excitations[0].get_magnetizing_current()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA,
            "The operating point carries no magnetizing current to read its peak from");
    }
    const auto magnetizingCurrent = excitations[0].get_magnetizing_current().value();
    const auto waveform = magnetizingCurrent.get_waveform();
    if (waveform && !waveform->get_data().empty()) {
        const auto data = waveform->get_data();
        return std::max(std::abs(*std::max_element(data.begin(), data.end())), std::abs(*std::min_element(data.begin(), data.end())));
    }
    if (magnetizingCurrent.get_processed() && magnetizingCurrent.get_processed()->get_peak()) {
        return std::abs(magnetizingCurrent.get_processed()->get_peak().value());
    }
    throw InvalidInputException(ErrorCode::MISSING_DATA,
        "The magnetizing current of the operating point has neither a waveform nor a processed peak");
}

// The peak magnetizing current this core and coil carry at an operating point: the one the
// magnetizing-inductance model stamps on it (a copy; the caller's operating point is not changed).
inline double model_magnetizing_current_peak(const Core& core, const Coil& coil, const Inputs& inputs, OperatingPoint operatingPoint) {
    MagnetizingInductance magnetizingInductanceModel;
    magnetizingInductanceModel.set_topology(magnetizing_model_topology(inputs));
    magnetizingInductanceModel.calculate_inductance_and_magnetic_flux_density(core, coil, &operatingPoint);
    return stamped_magnetizing_current_peak(operatingPoint);
}

// Single source of truth for "is this operating point a common-mode choke's",
// shared by the saturation FILTER (seed stage) and the MagneticAdviser's final
// saturation GATE so the two cannot disagree (ABT #1369: the gate lacked the
// exemption and dropped every wound CMC the filter had let through).
//
// A CMC has every winding on the line side, but its
// line current is DIFFERENTIAL: the flux of the line/neutral windings cancels
// in the core, so neither the raw line current nor the ripple swing that
// get_common_mode_choke_magnetizing_current models as "the CM current" is a
// saturation criterion — MAS carries no CM-noise-current requirement to gate
// against, and the ripple estimate puts multi-tesla fictional B on every EMI
// toroid. The impedance filter is the authoritative gate for a CMC. CMC-ONLY,
// deliberately: a differential-mode choke also has 2-4 matching-current
// windings (can_be_common_mode_choke is true for it too) but carries the full
// line current additively and genuinely saturates, so it keeps both gates.
// Distinguish by the declared role (declared_magnetic_role: the topology, or the
// web wizard's subApplication — the wizards do not set designRequirements.topology), then confirm current
// symmetry with can_be_common_mode_choke so a mis-tagged asymmetric design is
// still checked.
inline bool is_tagged_common_mode_choke(const Inputs& inputs, const OperatingPoint& operatingPoint) {
    auto declared = declared_magnetic_role(inputs);
    bool taggedCommonMode = declared && declared->role == MagneticRole::COMMON_MODE_CHOKE;
    return taggedCommonMode && Inputs::can_be_common_mode_choke(operatingPoint);
}

// Single source of truth for the saturation-safe "maximum allowed peak flux
// density" used as the SIZING target by the CoreAdviser gapping code
// (gap/turn selection, refinement, golden-section / binary-search) AND by the
// CoreAdviserDataset transformer turn seeder. Conservative ceiling: the smaller
// of the Maniktala saturation-margin limit (Bsat/margin) and the
// maximumProportion flux cap, so the gapping and the seeder size to exactly the
// SAME limit and cannot drift (previously the seeder used Bsat/margin = 0.833
// while the gapping capped at min(1/margin, maxProportion) = 0.7). `rawBsat`
// must be the proportion=false saturation flux at the operating temperature.
// The saturation filter independently rejects Bpeak·margin > Bsat.
inline double maximum_allowed_magnetic_flux_density(double rawBsat) {
    return rawBsat * std::min(1.0 / settings.get_core_adviser_saturation_margin(),
                              defaults.maximumProportionMagneticFluxDensitySaturation);
}

// Single source of truth for the temperature at which saturation (B_sat / I_sat)
// is EVALUATED across every CoreAdviser/MagneticAdviser saturation site (filter
// gate, transformer turn seeder, inductor turn seeder, fast-path final gate,
// loss-sweep cap, gap sizing). Ferrite B_sat falls with temperature, so the hot
// junction corner — not the 25 C ambient — is the real operating point and the
// one downstream derating checks use (Heaviside's isat sweep evaluates at 100 C).
// Returns max(operatingTemperature, derating floor) so a spec hotter than the
// floor is never made cooler. This is a TEMPERATURE derating, orthogonal to the
// multiplicative saturation margin — the two are deliberately NOT stacked with
// the 0.7 flux-proportion factor (ABT #13): the inductor gate uses RAW B_sat at
// this temperature, so its total derating is exactly (hot temperature × margin).
inline double saturation_derating_temperature(double operatingTemperature) {
    return std::max(operatingTemperature,
                    settings.get_core_adviser_saturation_derating_temperature());
}

// Phase 3 (F6): PQI / UI shapes use integrated windings whose layout is
// not produced by fast_wind(), so the loss / impedance filters take a
// per-shape policy branch. Centralised here so every site uses the same
// predicate.
inline bool is_pqi_or_ui_shape(const std::string& shapeName) {
    return shapeName.rfind("PQI", 0) == 0 || shapeName.rfind("UI ", 0) == 0;
}

// Phase 3 (F1): helpers shared between MagneticFilterCoreAndDcLosses and
// MagneticFilterCoreDcAndSkinLosses. These two classes are near-clones
// (same members, same constructors, same shaped evaluate_magnetic with
// minor sweep-step + skin-effect differences). Pulling the verbatim
// chunks out lets the two classes focus on what genuinely differs.

inline std::map<std::string, std::string> default_loss_filter_models() {
    std::map<std::string, std::string> models;
    models["gapReluctance"] = to_string(defaults.reluctanceModelDefault);
    models["coreLosses"] = to_string(defaults.coreLossesModelDefault);
    models["coreTemperature"] = to_string(defaults.coreTemperatureModelDefault);
    return models;
}

// ABT #825: the area filter and the loss-budget threshold both need the time-averaged |v*i| of
// the primary excitation, so they both need a voltage AND a current. Each dereferenced both
// optionals unconditionally, so an operating point carrying only a current -- an ordinary way to
// describe an inductor's winding, and what several catalogue callers actually send -- died on an
// empty optional. The message was "bad optional access": no operating point, no winding, no field,
// no hint that a voltage was what was missing. Core sizing genuinely cannot proceed without one,
// so this still fails -- it just says what it needed.
inline void require_voltage_and_current_for_power(const OperatingPointExcitation& excitation,
                                                  size_t operatingPointIndex,
                                                  const OperatingPoint& operatingPoint,
                                                  const std::string& consumer) {
    std::string label = "operating point " + std::to_string(operatingPointIndex) +
        (operatingPoint.get_name() ? " ('" + operatingPoint.get_name().value() + "')" : "");

    std::string missing;
    if (!excitation.get_voltage() || !excitation.get_voltage()->get_waveform()) {
        missing = "a voltage";
    }
    if (!excitation.get_current() || !excitation.get_current()->get_waveform()) {
        missing = missing.empty() ? "a current" : missing + " and a current";
    }
    if (missing.empty()) {
        return;
    }

    throw std::invalid_argument(
        label + ", primary winding: " + consumer + " needs the time-averaged |v*i| of the"
        " excitation, and this one has no " + missing + " waveform. Give the primary winding both a"
        " voltage and a current excitation.");
}


// Compute the max-over-operating-points of the time-averaged |v·i|
// "power mean" used to scale the loss-budget threshold. If any
// operating point has a waveform longer than 2× the configured sampled
// width, switch the core-losses model to Steinmetz (mutates `models`).
inline double compute_maximum_power_mean_and_maybe_force_steinmetz(
    Inputs& inputs, std::map<std::string, std::string>& models)
{
    bool largeWaveform = false;
    // Phase 6 (perf): cache operating-points by const-ref to avoid OperatingPoint deep copies.
    const auto& operatingPoints = inputs.get_operating_points();
    std::vector<double> powerMeans(operatingPoints.size(), 0);
    for (size_t opi = 0; opi < operatingPoints.size(); ++opi) {
        auto excitation = Inputs::get_primary_excitation(operatingPoints[opi]);
        require_voltage_and_current_for_power(excitation, opi, operatingPoints[opi],
                                             "the core-losses budget threshold");
        auto voltageWaveform = excitation.get_voltage().value().get_waveform().value();
        auto currentWaveform = excitation.get_current().value().get_waveform().value();
        double frequency = excitation.get_frequency();

        if (voltageWaveform.get_data().size() != currentWaveform.get_data().size()) {
            voltageWaveform = Inputs::calculate_sampled_waveform(voltageWaveform, frequency, std::max(voltageWaveform.get_data().size(), currentWaveform.get_data().size()));
            currentWaveform = Inputs::calculate_sampled_waveform(currentWaveform, frequency, std::max(voltageWaveform.get_data().size(), currentWaveform.get_data().size()));
        }
        std::vector<double> voltageWaveformData = voltageWaveform.get_data();
        std::vector<double> currentWaveformData = currentWaveform.get_data();
        if (currentWaveformData.size() > settings.get_inputs_number_points_sampled_waveforms() * 2) {
            largeWaveform = true;
        }
        for (size_t i = 0; i < voltageWaveformData.size(); ++i) {
            powerMeans[opi] += std::fabs(voltageWaveformData[i] * currentWaveformData[i]);
        }
        powerMeans[opi] /= voltageWaveformData.size();
    }
    if (largeWaveform) {
        models["coreLosses"] = to_string(CoreLossesModels::STEINMETZ);
    }
    return *std::max_element(powerMeans.begin(), powerMeans.end());
}

// Attach a quick bobbin to the magnetic and validate its window width.
// No-op for PQI / UI shapes (integrated-winding designs that don't take
// a separate bobbin).
inline void prepare_bobbin_for_non_pqi(Magnetic* magnetic, const std::string& shapeName) {
    if (is_pqi_or_ui_shape(shapeName)) return;
    auto bobbin = Bobbin::create_quick_bobbin(magnetic->get_core());
    magnetic->get_mutable_coil().set_bobbin(bobbin);
    auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
    if (windingWindows[0].get_width()) {
        if ((windingWindows[0].get_width().value() < 0) || (windingWindows[0].get_width().value() > 1)) {
            throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened in bobbins 1:   windingWindows[0].get_width(): " + std::to_string(static_cast<int>(bool(windingWindows[0].get_width()))) +
                                     " windingWindows[0].get_width().value(): " + std::to_string(windingWindows[0].get_width().value()) +
                                     " shapeName: " + shapeName);
        }
    }
}

// Pick Steinmetz vs Proprietary based on what the core advertises, then
// compute core losses. Returns {ok, output, value}:
//   ok=false  → tiny negative from Proprietary interpolation noise;
//              the caller should explicitly reject this candidate
//              (Phase 1 fix: was previously a silent `break`).
//   ok=true   → value is the non-negative loss; output is the full result.
// A genuinely negative Steinmetz result still throws via the caller's
// own `coreLosses < 0` check.
struct CoreLossesPick {
    bool ok;
    CoreLossesOutput output;
    double value;
};
inline CoreLossesPick compute_core_losses_with_negative_guard(
    const Core& core, OperatingPointExcitation& excitation, double temperature,
    const std::shared_ptr<CoreLossesModel>& steinmetz,
    const std::shared_ptr<CoreLossesModel>& proprietary)
{
    auto methods = core.get_available_core_losses_methods();
    CoreLossesPick pick{};
    if (std::find(methods.begin(), methods.end(), VolumetricCoreLossesMethodType::STEINMETZ) != methods.end()) {
        pick.output = steinmetz->get_core_losses(core, excitation, temperature);
        pick.value = pick.output.get_core_losses();
        pick.ok = true;
    } else {
        // The proprietary model only has formulas for a few manufacturers
        // (Magnetics/Micrometals/Poco/TDG). Materials characterised solely by a
        // loss factor or complex permeability — e.g. Fair-Rite suppression
        // ferrites (75/76/…) — have neither Steinmetz nor a proprietary formula,
        // so calling the proprietary model would throw MATERIAL_DATA_MISSING.
        // Treat such a candidate as loss-unscorable (ok=false) instead: the caller
        // rejects it, and the suppression pipeline's pass-through keeps it as a
        // candidate (losses aren't its ranking criterion).
        bool hasProprietaryFormula =
            std::find(methods.begin(), methods.end(), VolumetricCoreLossesMethodType::MAGNETICS) != methods.end() ||
            std::find(methods.begin(), methods.end(), VolumetricCoreLossesMethodType::MICROMETALS) != methods.end() ||
            std::find(methods.begin(), methods.end(), VolumetricCoreLossesMethodType::POCO) != methods.end() ||
            std::find(methods.begin(), methods.end(), VolumetricCoreLossesMethodType::TDG) != methods.end();
        if (!hasProprietaryFormula) {
            pick.ok = false;
            return pick;
        }
        pick.output = proprietary->get_core_losses(core, excitation, temperature);
        pick.value = pick.output.get_core_losses();
        pick.ok = pick.value >= 0;
    }
    return pick;
}

// Final scoring/emission block: compute the mean total loss across all
// operating points, compare against the maximum power mean budget,
// optionally emit per-OP outputs.
inline std::pair<bool, double> finalize_losses_scoring(
    const std::vector<double>& totalLossesPerOP,
    const std::vector<CoreLossesOutput>& coreLossesPerOP,
    const std::vector<WindingLossesOutput>& windingLossesPerOP,
    Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs,
    double maximumPowerMean)
{
    if (totalLossesPerOP.size() < inputs->get_operating_points().size()) {
        return {false, 0.0};
    }
    double meanTotalLosses = 0;
    for (size_t opi = 0; opi < inputs->get_operating_points().size(); ++opi) {
        meanTotalLosses += totalLossesPerOP[opi];
    }
    if (!std::isfinite(meanTotalLosses)) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happend in core losses calculation for magnetic: " + magnetic->get_manufacturer_info().value().get_reference().value());
    }
    meanTotalLosses /= inputs->get_operating_points().size();
    bool valid = meanTotalLosses < maximumPowerMean * defaults.coreAdviserMaximumPercentagePowerCoreLosses / defaults.coreAdviserThresholdValidity;

    if (outputs != nullptr) {
        for (size_t opi = 0; opi < inputs->get_operating_points().size(); ++opi) {
            while (outputs->size() < opi + 1) outputs->push_back(Outputs());
            (*outputs)[opi].set_core_losses(coreLossesPerOP[opi]);
            (*outputs)[opi].set_winding_losses(windingLossesPerOP[opi]);
        }
    }
    return {valid, meanTotalLosses};
}

// Phase 3 (F5): used by the LossesTimes*/VolumeTimes* combinator filters
// to fetch a previously-computed scoring out of the global scoring cache,
// or fall back to running the sub-filter inline if no cache hit. Replaces
// the seven-line "previousX = get_scoring(...); if (previousX) {...} else
// {compute}" block that appeared verbatim 7 times.
inline double cached_or_compute_scoring(Magnetic* magnetic, Inputs* inputs,
                                        std::vector<Outputs>* outputs,
                                        MagneticFilters cacheKey,
                                        MagneticFilter& fallbackFilter) {
    auto cached = get_scoring(magnetic->get_reference(), cacheKey);
    if (cached) return cached.value();
    return fallbackFilter.evaluate_magnetic(magnetic, inputs, outputs).second;
}

} // namespace OpenMagnetics
