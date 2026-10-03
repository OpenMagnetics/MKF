#include "advisers/MagneticFilter.h"
#include "advisers/MagneticFilterInternal.h"
#include "physical_models/CoreLosses.h"
#include "physical_models/LeakageInductance.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/Temperature.h"
#include "support/Exceptions.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace OpenMagnetics {

std::pair<bool, double> MagnetomotiveForce::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto coil = magnetic->get_coil();
    double maximumMagnetomotiveForce = 0;
    if (inputs->get_operating_points().size() > 0 && magnetic->get_mutable_coil().get_functional_description().size() != inputs->get_operating_points()[0].get_excitations_per_winding().size()) {
        return {false, 0.0};        
    }

    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
        std::vector<double> currentRmsPerParallelPerWinding;
        for (size_t windingIndex = 0; windingIndex < magnetic->get_mutable_coil().get_functional_description().size(); ++windingIndex) {
            auto excitation = inputs->get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex];
            if (!excitation.get_current()) {
                throw InvalidInputException("Current is missing in excitation");
            }
            if (!excitation.get_current()->get_processed()) {
                throw InvalidInputException("Current is not processed");
            }
            if (!excitation.get_current()->get_processed()->get_rms()) {
                throw InvalidInputException("Current RMS is not processed");
            }
            auto currentRms = excitation.get_current()->get_processed()->get_rms().value();
            currentRmsPerParallelPerWinding.push_back(currentRms / coil.get_functional_description()[windingIndex].get_number_parallels());
            if (!coil.get_layers_description()) {
                throw CoilNotProcessedException("Coil not wound");
            }
        }
        std::vector<double> magnetomotiveForcePerLayer;
        magnetomotiveForcePerLayer.push_back(0);
        auto layers = coil.get_layers_description().value();
        for (auto layer : layers) {
            double magnetomotiveForceThisLayer = magnetomotiveForcePerLayer.back();
            if (layer.get_type() == ElectricalType::CONDUCTION) {

                auto windingIndex = coil.get_winding_index_by_name(layer.get_partial_windings()[0].get_winding());
                auto numberTurns = coil.get_functional_description()[windingIndex].get_number_turns();  
                auto numberPhysicalTurnsInLayer = 0;
                for (auto parallelProportion : layer.get_partial_windings()[0].get_parallels_proportion()) {
                    numberPhysicalTurnsInLayer += round(numberTurns * parallelProportion);
                }
                numberPhysicalTurnsInLayer *= layer.get_partial_windings().size();
                if (coil.get_functional_description()[windingIndex].get_isolation_side() == IsolationSide::PRIMARY) {
                    magnetomotiveForceThisLayer += numberPhysicalTurnsInLayer * currentRmsPerParallelPerWinding[windingIndex];
                }
                else {
                    magnetomotiveForceThisLayer -= numberPhysicalTurnsInLayer * currentRmsPerParallelPerWinding[windingIndex];
                }
            }
            // Phase 1 fix: previously push_back happened only in the
            // INSULATION (else) branch, so the CONDUCTION branch's update
            // to magnetomotiveForceThisLayer was dropped (dead local) and
            // the per-layer MMF trace was always zero for conductors.
            magnetomotiveForcePerLayer.push_back(magnetomotiveForceThisLayer);
        }

        double maximumMagnetomotiveForceThisOperatingPoint = *max_element(magnetomotiveForcePerLayer.begin(), magnetomotiveForcePerLayer.end());
        double minimumMagnetomotiveForceThisOperatingPoint = *min_element(magnetomotiveForcePerLayer.begin(), magnetomotiveForcePerLayer.end());
        maximumMagnetomotiveForceThisOperatingPoint = std::max(fabs(maximumMagnetomotiveForceThisOperatingPoint), fabs(minimumMagnetomotiveForceThisOperatingPoint));
        maximumMagnetomotiveForce = std::max(maximumMagnetomotiveForce, maximumMagnetomotiveForceThisOperatingPoint);

    }
    return {true, maximumMagnetomotiveForce};
}

std::pair<bool, double> MagneticFilterLeakageInductance::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    switch (_mode) {
        case LeakageInductanceFilterMode::MINIMIZE_LEAKAGE_RATIO:
            return evaluate_minimize_leakage_ratio(magnetic, inputs, outputs);
        case LeakageInductanceFilterMode::TARGET:
            return evaluate_target(magnetic, inputs, outputs);
    }
    throw InvalidInputException("MagneticFilterLeakageInductance: unknown mode");
}

std::pair<bool, double> MagneticFilterLeakageInductance::evaluate_target(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    // Integrated-leakage designs (LLC resonant inductance in the transformer, designed flyback
    // leakage, shunted transformers): the leakage is a target band, not something to minimise.
    auto requirements = inputs->get_design_requirements().get_leakage_inductance();
    if (!requirements || requirements->empty()) {
        throw InvalidInputException("MagneticFilterLeakageInductance TARGET mode needs designRequirements.leakageInductance");
    }
    size_t numberWindings = magnetic->get_coil().get_functional_description().size();
    if (numberWindings < 2 || requirements->size() > numberWindings - 1) {
        throw InvalidInputException("MagneticFilterLeakageInductance TARGET mode: " + std::to_string(requirements->size()) +
                                    " leakage requirements for " + std::to_string(numberWindings) + " windings");
    }
    if (inputs->get_operating_points().empty()) {
        throw InvalidInputException("MagneticFilterLeakageInductance TARGET mode needs an operating point for the frequency");
    }
    double frequency = inputs->get_operating_points()[0].get_excitations_per_winding()[0].get_frequency();

    LeakageInductance leakageModel;
    bool valid = true;
    double scoring = 0;
    std::vector<DimensionWithTolerance> leakagePerWinding;
    std::string methodUsed;
    for (size_t requirementIndex = 0; requirementIndex < requirements->size(); ++requirementIndex) {
        auto requirement = requirements.value()[requirementIndex];
        auto leakageOutput = leakageModel.calculate_leakage_inductance(*magnetic, frequency, 0, requirementIndex + 1);
        methodUsed = leakageOutput.get_method_used();
        double leakageInductance = resolve_dimensional_values(leakageOutput.get_leakage_inductance_per_winding()[0]);
        leakagePerWinding.push_back(leakageOutput.get_leakage_inductance_per_winding()[0]);

        // Relative distance to the band: zero inside [minimum, maximum] (a missing bound is open),
        // measured from the nearest bound outside it; with only a nominal, from the nominal.
        double reference = resolve_dimensional_values(requirement);
        if (!(reference > 0)) {
            throw InvalidInputException("MagneticFilterLeakageInductance TARGET mode: leakage requirement " + std::to_string(requirementIndex) + " is not positive");
        }
        double distance = 0;
        if (requirement.get_minimum() || requirement.get_maximum()) {
            if (requirement.get_minimum() && leakageInductance < requirement.get_minimum().value()) {
                distance = requirement.get_minimum().value() - leakageInductance;
            }
            else if (requirement.get_maximum() && leakageInductance > requirement.get_maximum().value()) {
                distance = leakageInductance - requirement.get_maximum().value();
            }
        }
        else {
            distance = std::fabs(leakageInductance - reference);
        }
        scoring += distance / reference;
        if (!check_requirement(requirement, leakageInductance)) {
            valid = false;
        }
    }

    if (outputs != nullptr) {
        for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
            while (outputs->size() < operatingPointIndex + 1) {
                outputs->push_back(Outputs());
            }
            InductanceOutput inductanceOutput;
            if ((*outputs)[operatingPointIndex].get_inductance()) {
                inductanceOutput = *(*outputs)[operatingPointIndex].get_inductance();
            }
            LeakageInductanceOutput leakageOutput;
            leakageOutput.set_method_used(methodUsed);
            leakageOutput.set_origin(ResultOrigin::SIMULATION);
            leakageOutput.set_leakage_inductance_per_winding(leakagePerWinding);
            inductanceOutput.set_leakage_inductance(leakageOutput);
            (*outputs)[operatingPointIndex].set_inductance(inductanceOutput);
        }
    }
    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterLeakageInductance::evaluate_minimize_leakage_ratio(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    // Leakage inductance filter for CMC optimization
    // For Common Mode Chokes, we want to minimize leakage inductance to maximize coupling coefficient
    // Coupling coefficient k = 1 - (Lk / Lm), where Lk is leakage and Lm is magnetizing inductance
    // Lower leakage means higher coupling, which means better common mode rejection
    
    if (magnetic->get_coil().get_functional_description().size() < 2) {
        // Single winding - no leakage calculation possible
        return {true, 0.0};
    }

    // Use the first operating point frequency for leakage calculation
    double frequency = defaults.measurementFrequency;
    if (inputs->get_operating_points().size() > 0) {
        frequency = inputs->get_operating_points()[0].get_excitations_per_winding()[0].get_frequency();
    }

    // Phase 1 fix: previously this method swallowed every exception and
    // returned {false, DBL_MAX} as an in-band sentinel, and used the same
    // sentinel when magnetizingInductance was zero. Both hid real bugs
    // (the DBL_MAX value also poisoned any downstream normalization).
    // Let exceptions propagate and throw explicitly on zero Lm.
    LeakageInductance leakageModel;
    auto leakageOutput = leakageModel.calculate_leakage_inductance(*magnetic, frequency, 0, 1);
    double leakageInductance = resolve_dimensional_values(leakageOutput.get_leakage_inductance_per_winding()[0]);

    // Get magnetizing inductance for normalization
    OpenMagnetics::MagnetizingInductance magnetizingInductanceModel("ZHANG");
    OperatingPoint* operatingPoint = nullptr;
    if (inputs->get_operating_points().size() > 0) {
        operatingPoint = &inputs->get_mutable_operating_points()[0];
    }
    auto magnetizingOutput = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(
        magnetic->get_mutable_core(),
        magnetic->get_mutable_coil(),
        operatingPoint
    );
    double magnetizingInductance = resolve_dimensional_values(magnetizingOutput.get_magnetizing_inductance());

    if (magnetizingInductance <= 0) {
        throw InvalidInputException(
            "MagneticFilterLeakageInductance: magnetizing inductance is non-positive ("
            + std::to_string(magnetizingInductance) + " H), cannot compute leakage ratio");
    }

    // Score is leakage ratio Lk/Lm - lower is better for CMC.
    // A perfect CMC has k ≈ 1, meaning Lk/Lm ≈ 0.
    double scoring = leakageInductance / magnetizingInductance;

    if (outputs != nullptr && inputs->get_operating_points().size() > 0) {
        for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
            while (outputs->size() < operatingPointIndex + 1) {
                outputs->push_back(Outputs());
            }
            InductanceOutput inductanceOutput;
            if ((*outputs)[operatingPointIndex].get_inductance()) {
                inductanceOutput = *(*outputs)[operatingPointIndex].get_inductance();
            }
            inductanceOutput.set_leakage_inductance(leakageOutput);
            (*outputs)[operatingPointIndex].set_inductance(inductanceOutput);
        }
    }

    return {true, scoring};
}

MagneticFilterTemperature::MagneticFilterTemperature(Inputs inputs, double maximumTemperature)
    : _maximumTemperature(maximumTemperature)
{
}

std::pair<bool, double> MagneticFilterTemperature::evaluate_magnetic(
    Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs)
{
    // Phase 1 fix: previously wrapped in `try {} catch (...) { return {true,
    // _maximumTemperature}; }` which silently treated EVERY failure as a
    // pass at the threshold. That hid model bugs (e.g. null _coreLossesModel)
    // and produced "always-acceptable" temperature scores. Let exceptions
    // propagate so callers can either handle them or fail loudly.
    // A design has to survive its hottest operating point, so each one gets its own
    // thermal solve and the gate reads the maximum. This used to solve ONCE with the
    // core losses AVERAGED over the operating points (and the last one's ambient), which
    // let a part through whenever a cool point diluted a hot one (ABT #1412).
    if (inputs->get_operating_points().empty()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA,
            "Temperature filter needs at least one operating point to compute core losses");
    }
    const auto& core = magnetic->get_core();
    double maximumTemperature = -std::numeric_limits<double>::max();

    const auto& coil = magnetic->get_coil();
    const std::string magneticRef = inductance_flux_cache_key(*magnetic);

    // The thermal solve is the most expensive per-candidate step of the core stage, and the
    // MagneticAdviser calls the core adviser again with a larger request whenever too few cores
    // wind: the same candidates came back and were solved again (a 10 kW transformer spent over
    // seven minutes re-solving cores it had already judged). The verdict depends only on the
    // candidate (reference, material, turns, gaps) and on the operating points, so it is memoised
    // per advise (clear_scoring) under a key holding both.
    std::string temperatureCacheKey = "temperature|" + magneticRef;
    for (const auto& op : inputs->get_operating_points()) {
        temperatureCacheKey += "|" + json(op.get_conditions()).dump();
        for (const auto& excitation : op.get_excitations_per_winding()) {
            std::ostringstream excitationKey;
            excitationKey.precision(17);
            excitationKey << " f " << excitation.get_frequency();
            if (excitation.get_current() && excitation.get_current()->get_processed()) {
                const auto processed = excitation.get_current()->get_processed().value();
                excitationKey << " Irms " << processed.get_rms().value_or(-1) << " Ipk " << processed.get_peak().value_or(-1);
            }
            if (excitation.get_voltage() && excitation.get_voltage()->get_processed()) {
                const auto processed = excitation.get_voltage()->get_processed().value();
                excitationKey << " Vrms " << processed.get_rms().value_or(-1) << " Vpk " << processed.get_peak().value_or(-1);
            }
            temperatureCacheKey += excitationKey.str();
        }
    }
    if (auto cachedMaximumTemperature = get_scoring(temperatureCacheKey, MagneticFilters::TEMPERATURE_RISE)) {
        return {cachedMaximumTemperature.value() <= _maximumTemperature, cachedMaximumTemperature.value()};
    }

    // The copper: every winding (a stand-in coil is completed from the turns ratios), laid out
    // by fast_wind() for its per-turn losses.
    SettingsGuard<bool> coilDelimitGuard(settings,
        &Settings::get_coil_delimit_and_compact,
        &Settings::set_coil_delimit_and_compact, false);
    const std::string shapeName = core.get_shape_name();
    const bool windable = !is_pqi_or_ui_shape(shapeName);
    // Whether this copper fits the window is the window copper capacity screen's verdict (run
    // before this filter). The stand-in's layout is fast_wind()'s, with section widths in
    // proportion to copper area: a winding of few, thick turns can get a section narrower than
    // its conductor and overlap the bobbin by a fraction of a millimetre, a layout artifact the
    // coil stage does not share. Such a turn gets no conduction path to the enclosure on that
    // face (logged by the network, ABT #1454), which leaves the estimate on the hot side,
    // instead of throwing the whole candidate out of the advise.
    SettingsGuard<bool> strictGeometryGuard(settings,
        &Settings::get_thermal_network_strict_geometry,
        &Settings::set_thermal_network_strict_geometry, false);
    Magnetic lossesMagnetic = with_every_winding(*magnetic, *inputs);
    // The thermal network gets one node per laid-out conductor, and the stand-in carries each
    // winding's current in parallel skin-depth strands (15 per turn for 50 A at 50 kHz): a
    // 10 kW transformer solved hundreds of strand nodes per candidate, ~2 s each. The network
    // is solved on the same coil with each winding's strands merged into one round conductor
    // of their copper area (same turns, same copper, laid out as one bundle), and every merged
    // turn carries the losses of the strands it replaces.
    Magnetic thermalMagnetic = lossesMagnetic;
    if (windable) {
        if (std::holds_alternative<std::string>(lossesMagnetic.get_coil().get_bobbin())) {
            prepare_bobbin_for_non_pqi(&lossesMagnetic, shapeName);
        }
        thermalMagnetic = lossesMagnetic;
        auto windings = thermalMagnetic.get_coil().get_functional_description();
        const double wireTemperature = saturation_derating_temperature(inputs->get_operating_points()[0].get_conditions().get_ambient_temperature());
        for (auto& winding : windings) {
            if (winding.get_number_parallels() > 1) {
                auto strand = winding.resolve_wire();
                const double copperArea = strand.calculate_conducting_area() * winding.get_number_parallels();
                winding.set_wire(Wire::get_wire_for_conducting_area(copperArea, wireTemperature, true));
                winding.set_number_parallels(1);
            }
        }
        thermalMagnetic.get_mutable_coil().set_functional_description(windings);
        for (auto* laidOut : {&lossesMagnetic, &thermalMagnetic}) {
            laidOut->get_mutable_coil().fast_wind();
            if (!laidOut->get_coil().get_turns_description()) {
                // As in the loss filter: a coil fast_wind() cannot lay out is rejected explicitly.
                logEntry("Temperature filter: the coil of " + magnetic->get_reference() +
                         " cannot be laid out, so its copper cannot be placed in the thermal network; rejected",
                         "MagneticFilterTemperature", 2);
                add_scoring(temperatureCacheKey, MagneticFilters::TEMPERATURE_RISE, std::numeric_limits<double>::max());
                return {false, std::numeric_limits<double>::max()};
            }
        }
    }
    // Per-turn losses of the strand coil summed onto the merged coil's turns: the k-th turn of a
    // winding takes the losses of the k-th turn of each of its strands.
    auto merge_losses_per_turn = [&](const WindingLossesOutput& strandLosses) {
        const auto strandTurns = lossesMagnetic.get_coil().get_turns_description().value();
        const auto mergedTurns = thermalMagnetic.get_coil().get_turns_description().value();
        auto strandLossesPerTurn = strandLosses.get_winding_losses_per_turn();
        if (!strandLossesPerTurn || strandLossesPerTurn->size() != strandTurns.size()) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT,
                "Temperature filter: the winding losses of " + magnetic->get_reference() + " do not give one entry per laid-out turn");
        }
        auto coil = thermalMagnetic.get_coil();
        std::map<std::pair<size_t, size_t>, double> lossesPerWindingTurn;
        std::map<std::pair<size_t, size_t>, size_t> turnCounter;
        for (size_t turnIndex = 0; turnIndex < strandTurns.size(); ++turnIndex) {
            const auto& element = (*strandLossesPerTurn)[turnIndex];
            double losses = element.get_ohmic_losses() ? element.get_ohmic_losses()->get_losses() : 0.0;
            for (const auto& lossesPerElement : {element.get_skin_effect_losses(), element.get_proximity_effect_losses()}) {
                if (lossesPerElement) {
                    for (double harmonicLosses : lossesPerElement->get_losses_per_harmonic()) {
                        losses += harmonicLosses;
                    }
                }
            }
            const size_t windingIndex = coil.get_winding_index_by_name(strandTurns[turnIndex].get_winding());
            const size_t parallelIndex = strandTurns[turnIndex].get_parallel();
            const size_t turnInWinding = turnCounter[{windingIndex, parallelIndex}]++;
            lossesPerWindingTurn[{windingIndex, turnInWinding}] += losses;
        }
        std::vector<WindingLossesPerElement> mergedLossesPerTurn;
        std::map<size_t, size_t> mergedTurnCounter;
        for (const auto& turn : mergedTurns) {
            const size_t windingIndex = coil.get_winding_index_by_name(turn.get_winding());
            const size_t turnInWinding = mergedTurnCounter[windingIndex]++;
            auto found = lossesPerWindingTurn.find({windingIndex, turnInWinding});
            if (found == lossesPerWindingTurn.end()) {
                throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT,
                    "Temperature filter: turn " + turn.get_name() + " of " + magnetic->get_reference() + " has no strand turn to take its losses from");
            }
            OhmicLosses ohmicLosses;  // the turn's DC and skin-effect losses together
            ohmicLosses.set_losses(found->second);
            ohmicLosses.set_origin(ResultOrigin::SIMULATION);
            WindingLossesPerElement element;
            element.set_name(turn.get_name());
            element.set_ohmic_losses(ohmicLosses);
            mergedLossesPerTurn.push_back(element);
        }
        if (mergedLossesPerTurn.size() != lossesPerWindingTurn.size()) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT,
                "Temperature filter: the merged coil of " + magnetic->get_reference() + " lays out " + std::to_string(mergedLossesPerTurn.size()) +
                " turns for " + std::to_string(lossesPerWindingTurn.size()) + " strand turns");
        }
        WindingLossesOutput merged = strandLosses;
        merged.set_winding_losses_per_turn(mergedLossesPerTurn);
        return merged;
    };
    size_t opIndex = 0;
    for (auto& op : inputs->get_operating_points()) {
        double ambientTemperature = op.get_conditions().get_ambient_temperature();
        auto excitation = op.get_excitations_per_winding()[0];

        // Phase 8 (perf): if a prior filter (typically SATURATION on the
        // inductor branch) already ran calculate_inductance_and_magnetic_
        // flux_density for this (magnetic, OP) with the default reluctance
        // model, reuse its peak B. Both filters default-construct
        // MagnetizingInductance and use it on the same fixed coil, so the
        // peak B matches by construction. On cache miss, fall through to
        // the iterative path.
        SignalDescriptor magneticFluxDensity;
        if (auto cached = get_cached_inductance_flux(magneticRef, opIndex)) {
            magneticFluxDensity = cached->magneticFluxDensity;
        } else {
            auto opCopy = op;
            auto [L, B] = _magnetizingInductance.calculate_inductance_and_magnetic_flux_density(core, coil, &opCopy);
            (void)L;
            magneticFluxDensity = B;
        }
        excitation.set_magnetic_flux_density(magneticFluxDensity);
        // Delegate model selection to the orchestrator: the old
        // Steinmetz-else-proprietary dispatch sent every material without a
        // Steinmetz method to the proprietary model, which throws for
        // materials (e.g. Fair-Rite loss-factor ones) with no proprietary data
        CoreLossesOutput cl = _coreLosses.calculate_core_losses(core, excitation, ambientTemperature);

        TemperatureConfig config;
        config.coreLosses = cl.get_core_losses();
        if (windable) {
            // Copper resistivity at the same hot corner the loss filters evaluate at.
            const double windingTemperature = saturation_derating_temperature(ambientTemperature);
            auto windingLosses = _windingOhmicLosses.calculate_ohmic_losses(lossesMagnetic.get_coil(), op, windingTemperature);
            windingLosses = _windingSkinEffectLosses.calculate_skin_effect_losses(lossesMagnetic.get_coil(), windingTemperature, windingLosses,
                                                                                  settings.get_harmonic_amplitude_threshold());
            if (!(windingLosses.get_winding_losses() >= 0) || !std::isfinite(windingLosses.get_winding_losses())) {
                throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT,
                    "Temperature filter: invalid winding losses for " + magnetic->get_reference());
            }
            config.coreOnly = false;
            config.windingLosses = windingLosses.get_winding_losses();
            config.windingLossesOutput = merge_losses_per_turn(windingLosses);
        }
        else {
            config.coreOnly = true;
        }
        config.ambientTemperature = ambientTemperature;
        config.plotSchematic = false;
        if (op.get_conditions().get_cooling()) config.masCooling = op.get_conditions().get_cooling();

        Temperature temp(windable ? thermalMagnetic : *magnetic, config);
        auto result = temp.calculateTemperatures();
        maximumTemperature = std::max(maximumTemperature, result.maximumTemperature);
        ++opIndex;
    }

    add_scoring(temperatureCacheKey, MagneticFilters::TEMPERATURE_RISE, maximumTemperature);
    return {maximumTemperature <= _maximumTemperature, maximumTemperature};
}

} // namespace OpenMagnetics
