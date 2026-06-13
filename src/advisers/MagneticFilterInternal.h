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

// Phase 5: promoted from a static function inside MagneticFilter.cpp so the
// extracted filter translation units (MagneticFilterPhysical.cpp etc.) can
// share the topology classification without duplicating the switch.
// Returns true for energy-storing topologies (inductors, flyback-derived);
// false for transformer-style topologies and for unknown / unspecified
// values (treated as transformer, the safer no-gap default).
inline bool is_energy_storing_topology(std::optional<Topologies> topology) {
    if (!topology.has_value()) {
        return false; // Unknown topology, caller should use other heuristics
    }

    switch (topology.value()) {
        // Energy-storing topologies (inductors, flyback-derived)
        // These store energy in the magnetic field each switching cycle
        case Topologies::FLYBACK_CONVERTER:
        case Topologies::BUCK_CONVERTER:
        case Topologies::BOOST_CONVERTER:
        case Topologies::ISOLATED_BUCK_BOOST_CONVERTER:
        // Isolated Buck (Fly-Buck) is a synchronous buck with a coupled-inductor
        // secondary harvesting auxiliary outputs. The primary side stores energy
        // in the gap each switching cycle (buck inductor behavior). Industry
        // convention (TI AN-2292, Coilcraft) classifies the magnetic as a
        // "Fly-Buck coupled inductor", not a transformer -> needs a gapped core.
        case Topologies::ISOLATED_BUCK_CONVERTER:
        case Topologies::POWER_FACTOR_CORRECTION:
        // DMC is a single-winding (or balanced two-winding) inductor on the
        // line carrying the full DC line current as bias. Routing it through
        // the inductor B-from-current path lets MagnetizingInductance derate
        // permeability via the material's DC-bias polynomial; the transformer
        // path (B from voltage) misses the DC bias entirely.
        case Topologies::DIFFERENTIAL_MODE_CHOKE:
        // CMC is a two-winding choke that stores energy in the common-mode
        // inductance (Lcm) and carries the full DC line current as bias.
        // Like DMC, it must use the B-from-current path so that DC bias
        // derating is applied. The transformer path (B from voltage) is
        // inappropriate because CMCs are not voltage-driven transformers.
        case Topologies::COMMON_MODE_CHOKE:
        // SEPIC primary inductor L1 carries the full input DC current and stores
        // energy in the gap each switching cycle (buck-like). Even in the coupled-
        // inductor variant the magnetic returned by process_design_requirements is
        // L1 (see MAS schema sepic.json) -- still energy-storing. Needs a gap.
        case Topologies::SEPIC_CONVERTER:
        // Zeta is the dual of SEPIC: L1 is again a magnetizing/energy-storing
        // inductor (TI SLYT411 / SLVA721). Same classification as SEPIC.
        case Topologies::ZETA_CONVERTER:
        // Cuk: the input inductor L1 stores energy in the gap; the magnetic
        // produced by the wizard is L1 (or the coupled L1-L2 inductor).
        case Topologies::CUK_CONVERTER:
        case Topologies::FOUR_SWITCH_BUCK_BOOST_CONVERTER:
            return true;

        // Transformer topologies (forward-derived)
        // These transfer energy directly without storing it
        case Topologies::SINGLE_SWITCH_FORWARD_CONVERTER:
        case Topologies::TWO_SWITCH_FORWARD_CONVERTER:
        case Topologies::ACTIVE_CLAMP_FORWARD_CONVERTER:
        case Topologies::PUSH_PULL_CONVERTER:
        case Topologies::PHASE_SHIFTED_FULL_BRIDGE_CONVERTER:
        case Topologies::PHASE_SHIFTED_HALF_BRIDGE_CONVERTER:
        case Topologies::ASYMMETRIC_HALF_BRIDGE_CONVERTER:
        case Topologies::DUAL_ACTIVE_BRIDGE_CONVERTER:
        case Topologies::LLC_RESONANT_CONVERTER:
        case Topologies::CLLC_RESONANT_CONVERTER:
        case Topologies::CURRENT_TRANSFORMER:
            return false;

        default:
            return false; // Unknown, treat as transformer (safer - no gap)
    }
}

// A transformer must couple windings across at least two isolation sides. When every
// winding sits on the SAME isolation side (e.g. a coupled inductor with all windings on
// PRIMARY, like the Weinberg L1 input choke), the magnetic stores energy like an inductor
// regardless of the converter topology, so it must NOT be classified as a transformer.
// Misclassifying it as a transformer computes B from voltage (missing the DC bias), which
// over-estimates B and makes the saturation filter reject every core. A single (or unknown)
// winding returns false here so the topology / inductance heuristic decides as before.
inline bool windings_on_single_isolation_side(const std::optional<std::vector<IsolationSide>>& isolationSides) {
    if (!isolationSides.has_value() || isolationSides->size() < 2) {
        return false;
    }
    for (size_t i = 1; i < isolationSides->size(); ++i) {
        if ((*isolationSides)[i] != (*isolationSides)[0]) {
            return false;
        }
    }
    return true;
}

// Single source of truth for "is this design an energy-storing inductor"
// (as opposed to a transformer). This is THE classification that decides
// whether saturation scales with turns: for an inductor at a fixed gap more
// turns LOWER the saturation current (isat = B_sat·N·A_e/L, L ∝ N²), so the
// saturation-aware turn/gap co-design, the loss-sweep saturation cap, and the
// final realism gate all apply; for a transformer B is voltage-driven and more
// turns LOWER B, so they must not. Every one of those sites — plus the
// saturation filter and the CoreAdviser turn seeder — must use THIS predicate,
// or they disagree on which candidates to saturation-check. Detection tiers
// (most reliable first): (1) all windings on one isolation side ⇒ inductor;
// (2) topology, if specified; (3) legacy heuristic — a minimum-only inductance
// spec ("at least L") is a transformer's magnetizing inductance, anything with
// a nominal or maximum is a specific energy-storage target ⇒ inductor.
inline bool is_inductor(const Inputs& inputs) {
    const auto& designRequirements = inputs.get_design_requirements();
    if (windings_on_single_isolation_side(designRequirements.get_isolation_sides())) {
        return true;
    }
    auto topology = designRequirements.get_topology();
    if (topology.has_value()) {
        return is_energy_storing_topology(topology);
    }
    const auto& magnetizingInductance = designRequirements.get_magnetizing_inductance();
    bool minimumOnly = magnetizingInductance.get_minimum().has_value()
                       && !magnetizingInductance.get_nominal().has_value()
                       && !magnetizingInductance.get_maximum().has_value();
    return !minimumOnly;
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
