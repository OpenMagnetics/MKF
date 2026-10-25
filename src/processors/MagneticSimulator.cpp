#include "processors/MagneticSimulator.h"
#include "physical_models/LeakageInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "physical_models/Impedance.h"
#include "physical_models/ComplexPermeability.h"
#include "physical_models/Temperature.h"
#include "physical_models/InitialPermeability.h"
#include "support/Settings.h"
#include <MAS.hpp>
#include "support/Exceptions.h"
#include "support/NumberFormat.h"

#include <algorithm>
#include <memory>
#include <unordered_map>
#include <cmath>
#include <complex>
#include <limits>
#include <sstream>
#include <string>


namespace OpenMagnetics {

Mas MagneticSimulator::simulate(Mas mas, bool fastMode){
    return simulate(mas.get_mutable_inputs(), mas.get_mutable_magnetic(), fastMode);
}
Mas MagneticSimulator::simulate(const Inputs& inputs, const Magnetic& magnetic, bool fastMode){
    set_topology(magnetizing_model_topology(inputs));
    Mas mas;
    std::vector<Outputs> outputs;
    std::vector<OperatingPoint> simulatedOperatingPoints;
    mas.set_inputs(inputs);

    // Every winding must carry its own excitation (ABT #814). Deeper models
    // catch SOME of these — winding losses throw "Missing current excitation
    // for winding 2" once a later winding is unexcited — but a two-winding
    // magnetic given a single excitation used to run all the way through and
    // report complete, plausible losses for a component whose second winding
    // the caller never described. Silent under-specification is worse than a
    // refusal: the numbers look right and cannot be checked.
    {
        size_t numberWindings = magnetic.get_coil().get_functional_description().size();
        const auto& operatingPoints = mas.get_inputs().get_operating_points();
        for (size_t index = 0; index < operatingPoints.size(); ++index) {
            size_t numberExcitations = operatingPoints[index].get_excitations_per_winding().size();
            if (numberExcitations != numberWindings) {
                throw InvalidInputException(
                    ErrorCode::MISSING_DATA,
                    "Operating point " + std::to_string(index) + " has " +
                    std::to_string(numberExcitations) + " excitations for a magnetic with " +
                    std::to_string(numberWindings) + " windings; every winding needs exactly one excitation");
            }
        }
    }

    for (auto& operatingPoint : mas.get_mutable_inputs().get_mutable_operating_points()){
        Outputs output;
        // The full network needs turn nodes; wind a local copy if the coil is bare (as
        // calculate_temperature does).
        Magnetic woundMagnetic = magnetic;
        if (!woundMagnetic.get_coil().get_turns_description()) {
            woundMagnetic.get_mutable_coil().wind();
        }
        // ABT #2052 (Alf, 2026-10-24): at its ambient temperature the magnetizing model cannot read this point's
        // flux density on the material's magnetisation curve: it lies past the curve's end, above B_sat
        // (SaturationException) or, with a curve that ends below B_sat, under it (MagnetisationCurveExhaustedException).
        // The core is driven far past saturation there. That refusal affects only this point: its flux-dependent
        // outputs -- the magnetizing inductance, the flux density and core losses, and the temperature they heat
        // the part to -- are left unset (no value stands in for them), its winding losses are those of its
        // currents at the ambient temperature (as a point with no thermal operating point publishes them, ABT
        // #1876), and every other point is simulated. ResultsAnalyzer::analyze_results reports CORE_SATURATES for it.
        // At a heated iterate of the loss solve the same refusal is that point's NoThermalOperatingPointException
        // (ABT #1942); only at ambient does it reach here.
        auto simulateSaturatedPoint = [&]() {
            Outputs saturatedOutput;
            saturatedOutput.set_winding_losses(calculate_winding_losses(operatingPoint, woundMagnetic, std::nullopt));
            outputs.push_back(saturatedOutput);
            simulatedOperatingPoints.push_back(operatingPoint);
        };
        InductanceOutput inductanceOutput;
        // The model stamps its iterate's magnetizing current on the operating point before it reads the curve; a point
        // it refuses keeps its operating point as given.
        const OperatingPoint givenOperatingPoint = operatingPoint;
        try {
            inductanceOutput.set_magnetizing_inductance(calculate_magnetizing_inductance(operatingPoint, magnetic));
        }
        catch (const SaturationException&) {
            operatingPoint = givenOperatingPoint;
            simulateSaturatedPoint();
            continue;
        }
        catch (const MagnetisationCurveExhaustedException&) {
            operatingPoint = givenOperatingPoint;
            simulateSaturatedPoint();
            continue;
        }
        if (!fastMode) {
            if (magnetic.get_coil().get_functional_description().size() > 1) {
                inductanceOutput.set_leakage_inductance(calculate_leakage_inductance(operatingPoint, magnetic));
            }
        }
        output.set_inductance(inductanceOutput);
        // ABT #838: core and winding losses together, each at the temperature the COMBINED
        // dissipation produces (calculate_core_losses iterates both on one network). Fast mode
        // runs the same solve: the published core temperature is the same quantity in both.
        // ABT #1876 (Alf, 2026-10-14): a point with no thermal operating point (the network did not converge on
        // its losses, or its winding runs away) affects only itself: calculate_core_losses returns its core and
        // winding losses evaluated at ambient, as before ABT #838, and leaves temperatureOutput EMPTY, so it gets
        // NO temperature output. ResultsAnalyzer::analyze_results names it (NO_THERMAL_OPERATING_POINT). Any
        // other error refuses the whole design.
        std::optional<WindingLossesOutput> windingLossesOutput;
        std::optional<TemperatureOutput> temperatureOutput;
        output.set_core_losses(calculate_core_losses(operatingPoint, woundMagnetic, &windingLossesOutput, &temperatureOutput,
                                                     WindingLossesScope::FULL, NoThermalOperatingPointPolicy::LOSSES_AT_AMBIENT));
        output.set_winding_losses(windingLossesOutput.value());
        if (!fastMode && temperatureOutput) {
            // Full-network hot-spot from exactly the losses published above -- the same config
            // TemperatureConfig::fromSimulatedOutput builds, so the exported MAS carries the
            // temperature the UI's temperature map shows (ABT #906).
            output.set_temperature(temperatureOutput.value());
        }

        outputs.push_back(output);
        simulatedOperatingPoints.push_back(operatingPoint);
    }
    mas.get_mutable_inputs().set_operating_points(simulatedOperatingPoints);
    mas.set_magnetic(magnetic);
    mas.set_outputs(outputs);
    return mas;
}

MagnetizingInductanceOutput MagneticSimulator::calculate_magnetizing_inductance(OperatingPoint& operatingPoint, Magnetic magnetic){
    return _magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic, &operatingPoint);
}

LeakageInductanceOutput MagneticSimulator::calculate_leakage_inductance(OperatingPoint& operatingPoint, Magnetic magnetic){
    double frequency = operatingPoint.get_excitations_per_winding()[0].get_frequency();
    return calculate_leakage_inductance(magnetic, frequency);
}

LeakageInductanceOutput MagneticSimulator::calculate_leakage_inductance(Magnetic magnetic, double frequency){
    // Same shape as the public calculate_leakage_inductance WASM/Python API: one entry
    // per winding, index = winding index, 0 at the source (primary) slot, every value
    // referred to the primary. The previous hand-rolled loop emitted a secondaries-only
    // (N-1) array into the same MAS field, so consumers could not tell which shape they
    // were reading (web bug reports 125/139: the UI showed the primary's 0 as "Llk").
    return OpenMagnetics::LeakageInductance().calculate_leakage_inductance_all_windings(magnetic, frequency);
}

namespace {
    // ABT #1886: a coil described by its sections and turns but not its layers (the web Magnetic Builder keeps an
    // advised coil's sections and turns and drops its layers) has its layers rebuilt before the losses and the thermal
    // network read them: the winding losses group the turns by layer, and the network takes the insulation between
    // turns from the layers. Coil::rebuild_layers_description_from_turns throws when they cannot match the turns.
    void rebuild_layers_if_missing(Magnetic& magnetic) {
        const auto& coil = magnetic.get_coil();
        if (coil.get_turns_description() && coil.get_sections_description() && !coil.get_layers_description()) {
            magnetic.get_mutable_coil().rebuild_layers_description_from_turns();
        }
    }
}

WindingLossesOutput MagneticSimulator::calculate_winding_losses(OperatingPoint& operatingPoint, Magnetic magnetic, std::optional<double> temperature){
    auto& settings = OpenMagnetics::Settings::GetInstance();
    // RAII: if calculate_losses throws, the manual restore below never ran and the
    // process kept mirroring=1 for every subsequent computation.
    SettingsGuard<int> mirroringGuard(settings, &Settings::get_magnetic_field_mirroring_dimension, &Settings::set_magnetic_field_mirroring_dimension, 1);

    double simulationTemperature;
    if (!temperature) {
        simulationTemperature = operatingPoint.get_conditions().get_ambient_temperature();
    }
    else {
        simulationTemperature = temperature.value();
    }
    rebuild_layers_if_missing(magnetic);
    WindingLosses windingLosses;
    // windingLosses.set_winding_losses_harmonic_amplitude_threshold(0.01);
    auto losses = windingLosses.calculate_losses(magnetic, operatingPoint, simulationTemperature);
    return losses;
}

TemperatureOutput MagneticSimulator::calculate_temperature(OperatingPoint& operatingPoint, Magnetic magnetic, const Outputs& output) {
    // The full thermal network needs turn nodes; wind a local copy if the coil is bare.
    if (!magnetic.get_coil().get_turns_description()) {
        magnetic.get_mutable_coil().wind();
    }
    rebuild_layers_if_missing(magnetic);
    auto config = TemperatureConfig::fromSimulatedOutput(operatingPoint, output);
    auto thermalResult = Temperature(magnetic, config).calculateTemperatures();

    TemperatureOutput temperatureOutput;
    temperatureOutput.set_initial_temperature(config.ambientTemperature);
    temperatureOutput.set_maximum_temperature(thermalResult.maximumTemperature);
    double totalLosses = config.coreLosses + config.windingLosses;
    if (totalLosses > 0) {
        temperatureOutput.set_bulk_thermal_resistance(
            (thermalResult.maximumTemperature - config.ambientTemperature) / totalLosses);
    }
    temperatureOutput.set_method_used("ThermalNetwork");
    temperatureOutput.set_origin(ResultOrigin::SIMULATION);
    return temperatureOutput;
}

namespace {
    struct LossSolveMemo {
        std::unordered_map<std::string, SolvedLosses> entries;
        size_t hits = 0;
        size_t misses = 0;
    };
    thread_local std::unique_ptr<LossSolveMemo> activeLossSolveMemo;
}

LossSolveMemoScope::LossSolveMemoScope() {
    if (!activeLossSolveMemo) {
        activeLossSolveMemo = std::make_unique<LossSolveMemo>();
        _owner = true;
    }
}

LossSolveMemoScope::~LossSolveMemoScope() {
    if (_owner) {
        activeLossSolveMemo.reset();
    }
}

LossSolveMemoScope::Statistics LossSolveMemoScope::get_statistics() {
    if (!activeLossSolveMemo) {
        return {};
    }
    return {activeLossSolveMemo->hits, activeLossSolveMemo->misses, activeLossSolveMemo->entries.size()};
}

std::string MagneticSimulator::loss_solve_key(const OperatingPoint& operatingPoint, const Magnetic& magnetic,
                                              WindingLossesScope windingLossesScope) const {
    json magneticJson;
    to_json(magneticJson, magnetic);
    json operatingPointJson;
    to_json(operatingPointJson, operatingPoint);
    json modelsJson;
    modelsJson["windingLossesScope"] = static_cast<int>(windingLossesScope);
    modelsJson["magnetizingInductance"] = _magnetizingInductanceModel.get_models();
    std::vector<int> coreLossesModels;
    for (auto model : _coreLossesModel.get_core_losses_model_names()) {
        coreLossesModels.push_back(static_cast<int>(model));
    }
    modelsJson["coreLosses"] = coreLossesModels;
    auto requestedModel = _coreLossesModel.get_core_losses_requested_model_name();
    modelsJson["coreLossesRequested"] = requestedModel ? static_cast<int>(requestedModel.value()) : -1;
    // The effective value: the fingerprint holds what the caller set, which an adviser overrides.
    modelsJson["lossTemperatureLoop"] = Settings::GetInstance().get_simulator_loss_temperature_loop();
    return magneticJson.dump() + "\n" + operatingPointJson.dump() + "\n" + modelsJson.dump() + "\n" +
           Settings::GetInstance().fingerprint();
}

CoreLossesOutput MagneticSimulator::calculate_core_losses(OperatingPoint& operatingPoint, Magnetic magnetic,
                                                          std::optional<WindingLossesOutput>* windingLossesOutput,
                                                          std::optional<TemperatureOutput>* temperatureOutput,
                                                          WindingLossesScope windingLossesScope,
                                                          NoThermalOperatingPointPolicy noThermalOperatingPointPolicy) {
    if (!operatingPoint.get_excitations_per_winding()[0].get_current()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Missing current in operating point");
    }
    if (windingLossesOutput != nullptr && !magnetic.get_coil().get_turns_description()) {
        throw CoilNotProcessedException("Winding losses asked of '" + magnetic.get_reference() +
                                        "', whose coil has no turns: wind it first");
    }
    rebuild_layers_if_missing(magnetic);

    _lastDcBiasPoint.reset();
    _lastLowestFixedPointSearchEvaluations = 0;
    SolvedLosses solved;
    try {
        if (activeLossSolveMemo) {
            std::string key = loss_solve_key(operatingPoint, magnetic, windingLossesScope);
            auto found = activeLossSolveMemo->entries.find(key);
            if (found != activeLossSolveMemo->entries.end()) {
                ++activeLossSolveMemo->hits;
                solved = found->second;
            }
            else {
                ++activeLossSolveMemo->misses;
                solved = solve_losses_at_operating_temperature(operatingPoint, magnetic, windingLossesScope);
                activeLossSolveMemo->entries.emplace(std::move(key), solved);
            }
        }
        else {
            solved = solve_losses_at_operating_temperature(operatingPoint, magnetic, windingLossesScope);
        }
    }
    catch (const NoThermalOperatingPointException& noThermalOperatingPoint) {
        if (noThermalOperatingPointPolicy == NoThermalOperatingPointPolicy::THROW) {
            throw;
        }
        // ABT #1876 (Alf, 2026-10-14): this point has no thermal operating point. Its losses are those evaluated
        // at ambient (the solve's first iterate) and it has NO temperature: temperatureOutput is left empty.
        const auto& lossesAtAmbient = noThermalOperatingPoint.get_losses_at_ambient();
        operatingPoint.get_mutable_excitations_per_winding()[0] = lossesAtAmbient.excitation;
        if (windingLossesOutput != nullptr) {
            *windingLossesOutput = lossesAtAmbient.windingLosses;
        }
        if (temperatureOutput != nullptr) {
            *temperatureOutput = std::nullopt;
        }
        return lossesAtAmbient.coreLosses;
    }

    // Publish exactly the set that was solved.
    _lastDcBiasPoint = solved.dcBiasPoint;
    _lastLowestFixedPointSearchEvaluations = solved.lowestFixedPointSearchEvaluations;
    operatingPoint.get_mutable_excitations_per_winding()[0] = solved.excitation;
    if (windingLossesOutput != nullptr) {
        *windingLossesOutput = solved.windingLosses;
    }
    if (temperatureOutput != nullptr) {
        *temperatureOutput = solved.temperature;
    }
    return solved.coreLosses;
}

Magnetic MagneticSimulator::make_paramagnetic_magnetic(const Magnetic& magnetic) {
    // The magnetic with its core material made paramagnetic, mu_r = 1 (ABT #1855, Alf 2026-10-14): what the
    // magnetizing model sees of a core held at its Curie temperature (see solve_losses_at_operating_temperature).
    // Only the permeability changes: geometry, gaps and the saturation data stay; the core losses are still
    // evaluated on the real material.
    if (std::holds_alternative<std::vector<MaterialElement>>(magnetic.get_core().get_functional_description().get_material())) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT,
            "Losses of '" + magnetic.get_reference() + "': the losses drive its core past the temperatures its material data "
            "covers, and a core of more than one grade cannot be held there as a paramagnetic core");
    }
    Magnetic paramagnetic = magnetic;
    auto coreMaterial = paramagnetic.get_core().resolve_material();
    PermeabilityPoint vacuumPermeability;
    vacuumPermeability.set_value(1);
    auto permeabilities = coreMaterial.get_permeability();
    permeabilities.set_initial(vacuumPermeability);
    permeabilities.set_amplitude(std::nullopt);
    permeabilities.set_complex(std::nullopt);
    permeabilities.set_incremental(std::nullopt);
    coreMaterial.set_permeability(permeabilities);
    // Its own name: caches keyed by a material's name must not mix it with the ferromagnetic record.
    coreMaterial.set_name(coreMaterial.get_name() + " (paramagnetic, mu_r = 1)");
    paramagnetic.get_mutable_core().get_mutable_functional_description().set_material(coreMaterial);
    return paramagnetic;
}

bool MagneticSimulator::is_core_held_at_data_ceiling(const CoreMaterial& coreMaterial, double coreLossesTemperature, double maximumTemperature) {
    return is_core_held_at_data_ceiling(get_ferromagnetic_data_ceiling(coreMaterial), coreLossesTemperature, maximumTemperature);
}

bool MagneticSimulator::is_core_held_at_data_ceiling(const std::optional<CoreDataCeiling>& ceiling, double coreLossesTemperature,
                                                     double maximumTemperature) {
    // Held: the losses drove the core to its ceiling and the solve evaluated them there, within the loop's
    // tolerance below it. Until ABT #2078 this also read the hottest point, because a core held at Curie as the
    // paramagnetic core it is there loses almost nothing and the network came back cold (a P 18/11 TM63 resonant
    // inductor at 0.37 T: 23 W at 25 C, held at 289 C, then 0 W and a 33 C hottest point; Henry, 2026-10-24).
    // The solve now publishes at least the held temperature as the hottest point, so it carries no sign of the
    // hold. A point whose core losses were evaluated within the tolerance of the data limit is at that limit
    // whatever the network gives, and is reported as such.
    (void) maximumTemperature;
    return ceiling && coreLossesTemperature > ceiling->temperature - kOperatingTemperatureTolerance &&
           coreLossesTemperature <= ceiling->temperature;
}

std::optional<CoreDataCeiling> MagneticSimulator::get_ferromagnetic_data_ceiling(const CoreMaterial& coreMaterial) {
    const auto ceiling = InitialPermeability::get_maximum_ferromagnetic_data_temperature(coreMaterial);
    if (!ceiling) {
        return std::nullopt;
    }
    const auto curieTemperature = coreMaterial.get_curie_temperature();
    const bool isCurie = curieTemperature && ceiling.value() == curieTemperature.value();
    return CoreDataCeiling{ceiling.value(), isCurie ? CoreDataCeilingSource::CURIE_TEMPERATURE
                                                    : CoreDataCeilingSource::PERMEABILITY_TABLE_END};
}

std::optional<CoreDataCeiling> MagneticSimulator::get_core_data_ceiling(const Core& core, const OperatingPointExcitation& excitation,
                                                                        double ambientTemperature, CoreLosses& coreLossesModel) {
    const auto coreMaterial = core.resolve_material();
    auto ceiling = get_ferromagnetic_data_ceiling(coreMaterial);
    auto keepLower = [&](std::optional<double> end, CoreDataCeilingSource source) {
        if (end && (!ceiling || end.value() < ceiling->temperature)) {
            ceiling = CoreDataCeiling{end.value(), source};
        }
    };
    // H_sat(T) needs no part of its own: past its data it goes linearly to zero at the Curie temperature
    // (Core::get_magnetic_field_strength_saturation, Alf 2026-10-17), so it is positive below the Curie part.
    keepLower(coreLossesModel.get_temperature_factor_data_end(core, excitation, ambientTemperature),
              CoreDataCeilingSource::STEINMETZ_TEMPERATURE_FACTOR_END);
    return ceiling;
}

SolvedLosses MagneticSimulator::solve_losses_at_operating_temperature(const OperatingPoint& operatingPoint,
                                                                                         Magnetic& magnetic,
                                                                                         WindingLossesScope windingLossesScope) {
    // ABT #838: a wound coil heats the core, and the core's temperature sets the winding's: both
    // losses are evaluated at the temperatures their COMBINED dissipation produces, on the full
    // network (turn nodes carry the winding power). Without turns there are no winding losses to
    // compute and no turn nodes to put them in (the core advisers' sentinel-coil candidates): the
    // core is then solved on the core-only network with its own losses.
    const bool withWinding = magnetic.get_coil().get_turns_description().has_value();
    // Settings::simulatorLossTemperatureLoop off (a user choice for simulate and the builder; advisers read it
    // on, see Settings::get_simulator_loss_temperature_loop): no loop. The first iterate below IS the one pass --
    // core and winding losses at the operating temperature, the network run once on them for the hottest
    // point -- and it is returned as it is, its temperature fed back into nothing (Alf, 2026-10-22).
    const bool lossTemperatureLoop = Settings::GetInstance().get_simulator_loss_temperature_loop();

    const double ambientTemperature = operatingPoint.get_conditions().get_ambient_temperature();
    double coreTemperature = ambientTemperature;
    double windingTemperature = ambientTemperature;
    std::string history;
    double coreRelaxation = 1.0;
    double windingRelaxation = 1.0;
    std::optional<double> previousCoreResidual;
    std::optional<double> previousWindingResidual;
    // The loop's ceiling (get_core_data_ceiling): the highest temperature where the material is
    // ferromagnetic, its mu_i(T) has data (the last table point when the table ends below Curie, else
    // Curie), and the Steinmetz ranges its losses evaluate have a positive ct(T) (ABT #1908). Above it the magnetizing model or the core-loss model refuses
    // (ABT #1739), so no iterate may land there. It is taken on the first iterate's excitation, once
    // its losses have been evaluated (so a material or frequency the losses cannot price throws there,
    // as it always did); the frequencies it reads do not change between iterates.
    std::optional<CoreDataCeiling> dataCeiling;
    bool dataCeilingTaken = false;
    // The winding losses' field step (H at every turn, and the current split between a winding's
    // parallels that feeds it, ABT #1314) is computed ONCE, at ambient, and every iterate runs only
    // the per-turn kernels at its own winding temperature (the architecture's two steps): the
    // conductors' resistance follows the iterate, the field and the split do not (decided by Alf
    // 2026-10-13). Without parallels H does not depend on the temperature at all; with parallels
    // the losses differ from a field redone at the solved temperature by about 1e-5 relative.
    WindingLosses windingLossesModel;
    std::optional<WindingLosses::WindingLossesField> windingLossesField;
    // Set once the losses drive the core to that ceiling (see below): from then on the core is held at
    // the last temperature its material's data covers and only the winding temperature is iterated.
    // When that ceiling is the Curie temperature the held core is magnetized as the paramagnetic core it
    // is above it (mu_r = 1, paramagneticMagnetic); when it is the end of a mu_i(T) table that stops
    // below Curie, the core is still ferromagnetic there and keeps its last tabulated permeability.
    bool coreHeldAtCeiling = false;
    std::optional<Magnetic> paramagneticMagnetic;
    // The previous iterate of the winding-only loop (core held), for the runaway check below.
    struct HeldWindingIterate {
        double windingTemperature;
        double nextWindingTemperature;
        double windingLosses;
    };
    std::optional<HeldWindingIterate> heldWindingIterate;
    // The first iterate's losses: core and winding at ambient (ABT #1876), what a point with no thermal
    // operating point publishes.
    std::optional<LossesAtAmbient> lossesAtAmbient;

    // One evaluation of the loop's map: the losses with the core at coreAt and the winding at windingAt, and the
    // network they drive. The iteration below and the search for its lowest stable fixed point (ABT #1840) both
    // read it, so they evaluate the same map.
    struct LoopEvaluation {
        MAS::SignalDescriptor magneticFluxDensity;
        std::optional<DcBiasPoint> dcBiasPoint;
        OperatingPointExcitation excitation;
        CoreLossesOutput coreLossesOutput;
        std::optional<WindingLossesOutput> windingLosses;
        TemperatureConfig config;
    };
    // ABT #1942 (Alf, 2026-10-19): the magnetizing model reads the current of a flux given by the voltage back from the
    // material's magnetisation curve at the core's temperature, and refuses where that curve cannot carry the flux: past
    // its end above B_sat (SaturationException: the core is driven far past saturation there), or past its end while
    // still under B_sat (MagnetisationCurveExhaustedException: near the Curie point H_sat(T) goes to zero and the curve
    // ends low). At a temperature the loss loop's iterates drove the core to, that point has no thermal operating point
    // the models can evaluate: the same answer as a loop that does not settle (ABT #1876, #1960), with its losses at
    // ambient. Before any heating (the first iterate, at ambient) the refusal is the design's own and propagates.
    auto throwNoThermalOperatingPoint = [&](double coreAt, const std::string& reason) {
        if (!lossesAtAmbient) {
            throw;
        }
        throw NoThermalOperatingPointException(ErrorCode::CALCULATION_INVALID_RESULT, "No thermal operating point: at a core "
            "temperature of " + std::to_string(coreAt) + " C, where the loss loop's iterates drove it, the magnetizing model "
            "cannot read this point's current on the material's curve: " + reason + " History:" + history,
            lossesAtAmbient.value());
    };
    auto evaluateLosses = [&](double coreAt, double windingAt) {
        LoopEvaluation evaluation;
        // Core losses at the core's temperature. The magnetizing model reads the core temperature
        // from the operating point's conditions (permeability vs T sets B for a current-driven
        // part), so it gets a COPY at T_core; the caller's operating point keeps its real ambient.
        OperatingPoint operatingPointAtCoreTemperature = operatingPoint;
        operatingPointAtCoreTemperature.get_mutable_conditions().set_ambient_temperature(coreAt);
        try {
            evaluation.magneticFluxDensity = _magnetizingInductanceModel.calculate_inductance_and_magnetic_flux_density(
                paramagneticMagnetic ? paramagneticMagnetic.value() : magnetic, &operatingPointAtCoreTemperature).second;
        }
        catch (const MagnetisationCurveExhaustedException& error) {
            throwNoThermalOperatingPoint(coreAt, error.message());
        }
        catch (const SaturationException& error) {
            throwNoThermalOperatingPoint(coreAt, error.message());
        }
        // A paramagnetic core has no knee for a bias to drive it past.
        if (!paramagneticMagnetic) {
            evaluation.dcBiasPoint = _magnetizingInductanceModel.get_last_dc_bias_point();
        }
        // Read the excitation after the inductance model: it stamps the magnetizing current it
        // derived (a CMC's common-mode current, a voltage's integral) into the operating point it
        // was given. A copy taken before it lacked that current whenever the operating point
        // arrived without one, and the loss-factor model then refused it ("Missing magnetizing current").
        evaluation.excitation = operatingPointAtCoreTemperature.get_excitations_per_winding()[0];
        evaluation.excitation.set_magnetic_flux_density(evaluation.magneticFluxDensity);
        evaluation.coreLossesOutput = _coreLossesModel.calculate_core_losses(magnetic.get_core(), evaluation.excitation, coreAt);
        evaluation.coreLossesOutput.set_temperature(coreAt);
        if (!dataCeilingTaken) {
            dataCeiling = get_core_data_ceiling(magnetic.get_core(), evaluation.excitation, ambientTemperature, _coreLossesModel);
            dataCeilingTaken = true;
        }

        // One network, every heat source that is modelled.
        if (withWinding) {
            switch (windingLossesScope) {
                case WindingLossesScope::FULL: {
                    // As calculate_winding_losses: the field mirrored in one dimension.
                    auto& settings = Settings::GetInstance();
                    SettingsGuard<int> mirroringGuard(settings, &Settings::get_magnetic_field_mirroring_dimension,
                                                      &Settings::set_magnetic_field_mirroring_dimension, 1);
                    if (!windingLossesField) {
                        windingLossesField = windingLossesModel.calculate_magnetic_field(magnetic, operatingPoint, ambientTemperature);
                    }
                    evaluation.windingLosses = windingLossesModel.calculate_losses(magnetic, windingLossesField.value(), windingAt);
                    break;
                }
                case WindingLossesScope::OHMIC_AND_SKIN_EFFECT:
                    evaluation.windingLosses = WindingSkinEffectLosses::calculate_skin_effect_losses(
                        magnetic.get_coil(), windingAt,
                        WindingOhmicLosses::calculate_ohmic_losses(magnetic.get_coil(), operatingPoint, windingAt),
                        0.5);
                    break;
                case WindingLossesScope::OHMIC:
                    evaluation.windingLosses = WindingOhmicLosses::calculate_ohmic_losses(magnetic.get_coil(), operatingPoint, windingAt);
                    break;
            }
            Outputs losses;
            losses.set_core_losses(evaluation.coreLossesOutput);
            losses.set_winding_losses(evaluation.windingLosses.value());
            evaluation.config = TemperatureConfig::fromSimulatedOutput(operatingPoint, losses);
        }
        else {
            evaluation.config = TemperatureConfig::fromMasOperatingConditions(operatingPoint.get_conditions());
            evaluation.config.plotSchematic = false;
            evaluation.config.coreLosses = evaluation.coreLossesOutput.get_core_losses();
            evaluation.config.coreOnly = true;
        }
        return evaluation;
    };
    // Whether the search below has run: once per solve, and the evaluations it made (SolvedLosses).
    bool lowestFixedPointSearched = false;
    size_t lowestFixedPointSearchEvaluations = 0;
    // h at ambient, from the first iterate: its core and winding losses are both evaluated at ambient, which is the
    // diagonal point the search would otherwise evaluate again.
    std::optional<std::pair<double, double>> residualAtAmbient;
    // The loop's map along the diagonal, h(T) = g(T) - T, where g(T) is the core temperature the network returns
    // for the core and winding losses evaluated at T, and the winding temperature it returns with it.
    auto coreResidualAt = [&](double temperature, size_t iteration) {
        ++lowestFixedPointSearchEvaluations;
        auto evaluation = evaluateLosses(temperature, temperature);
        ThermalResult thermalResult;
        try {
            thermalResult = Temperature(magnetic, evaluation.config).calculateTemperatures();
        }
        catch (const ThermalNetworkNotConvergedException& error) {
            throw NoThermalOperatingPointException(error.code(), "No thermal operating point: " + error.message() +
                " Loss loop at iterate " + std::to_string(iteration) + ", searching for its lowest stable fixed point (losses at " +
                format_significant(temperature) + " C). History:" + history, lossesAtAmbient.value());
        }
        if (!thermalResult.maximumCoreTemperature || (withWinding && !thermalResult.averageTurnTemperature)) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT,
                "Thermal network of '" + magnetic.get_reference() + "' has no core or turn node: no temperature to search the loss loop's fixed point with");
        }
        const double windingAt = withWinding ? thermalResult.averageTurnTemperature.value() : ambientTemperature;
        return std::make_pair(thermalResult.maximumCoreTemperature.value() - temperature, windingAt);
    };
    // The lowest T in (lowest, highest] where h falls through zero (h >= 0 below, h < 0 above: g' < 1 there, a stable
    // fixed point), and the winding temperature the network gives there; h(lowest) is atLowest. highest is the loop's
    // own iterate, within its tolerance of the data ceiling (the ceiling itself is where the models stop).
    // A coarse pass walks up the bracket in n equal steps, highest included, and stops at the first point with h < 0;
    // that sign change alone is then bisected to the loop's tolerance. n = ceil(sqrt(width / tolerance)): the spacing
    // is the geometric mean of the bracket and the tolerance, so a part with no fixed point (the common case here)
    // pays about the square root of a tolerance-spaced scan (ETD 29 N87, 25 C to its 210 C ceiling: 10 networks instead of the 63 of a scan refined to 2 K),
    // and one that has a fixed point pays its coarse steps plus log2(spacing / tolerance) bisections. highest is
    // evaluated because h can be negative there although the loop climbed to it: its iterates are not on the
    // diagonal (their winding losses lag), and a stable point can sit within the last step below the ceiling (a
    // 3C91 E 114 coil, 25 C: h < 0 from 209 C to its 220 C ceiling). What it cannot see is a stretch of h < 0 that
    // lies strictly between two coarse points: a stable and an unstable fixed point closer together than the
    // spacing, which happens only next to the load where the two merge and the stable point disappears. Nothing
    // (std::nullopt) when no coarse point has h < 0.
    auto findLowestStableCoreTemperature = [&](double lowest, double highest, std::pair<double, double> atLowest,
                                               size_t iteration) -> std::optional<std::pair<double, double>> {
        if (atLowest.first <= 0) {
            return std::make_pair(lowest, atLowest.second);
        }
        if (highest - lowest <= kOperatingTemperatureTolerance) {
            // No room above the lowest point for a fixed point the loop could tell from it.
            return std::nullopt;
        }
        const size_t steps = static_cast<size_t>(std::ceil(std::sqrt((highest - lowest) / kOperatingTemperatureTolerance)));
        const double spacing = (highest - lowest) / static_cast<double>(steps);
        double hotBelow = lowest;
        for (size_t index = 1; index <= steps; ++index) {
            const double temperature = index == steps ? highest : lowest + static_cast<double>(index) * spacing;
            auto atTemperature = coreResidualAt(temperature, iteration);
            if (atTemperature.first >= 0) {
                hotBelow = temperature;
                continue;
            }
            double coldAbove = temperature;
            while (coldAbove - hotBelow > kOperatingTemperatureTolerance) {
                const double middle = (hotBelow + coldAbove) / 2;
                auto atMiddle = coreResidualAt(middle, iteration);
                if (atMiddle.first < 0) {
                    coldAbove = middle;
                    atTemperature = atMiddle;
                }
                else {
                    hotBelow = middle;
                }
            }
            return std::make_pair(coldAbove, atTemperature.second);
        }
        return std::nullopt;
    };

    for (size_t iteration = 1; iteration <= kOperatingTemperatureMaximumIterations; ++iteration) {
        auto evaluation = evaluateLosses(coreTemperature, windingTemperature);
        auto& magneticFluxDensity = evaluation.magneticFluxDensity;
        const auto& dcBiasPoint = evaluation.dcBiasPoint;
        auto& excitation = evaluation.excitation;
        auto& coreLossesOutput = evaluation.coreLossesOutput;
        auto& windingLosses = evaluation.windingLosses;
        auto& config = evaluation.config;
        if (!lossesAtAmbient) {
            lossesAtAmbient = LossesAtAmbient{coreLossesOutput, excitation, windingLosses};
        }
        ThermalResult thermalResult;
        try {
            thermalResult = Temperature(magnetic, config).calculateTemperatures();
        }
        catch (const ThermalNetworkNotConvergedException& error) {
            if (!lossTemperatureLoop) {
                // One pass: there is no loop whose operating point could be missing. The network gives no
                // temperature for these losses, and the one-pass estimate cannot be made.
                throw ThermalNetworkNotConvergedException(error.message() + " One-pass temperature estimate of '" +
                    magnetic.get_reference() + "' (loss-temperature loop off): core and winding losses at the operating "
                    "temperature, " + std::to_string(coreTemperature) + " C.");
            }
            // No steady temperature exists on the network for this iterate's losses: the operating point has
            // no thermal operating point (ABT #1876, Alf 2026-10-14).
            throw NoThermalOperatingPointException(error.code(), "No thermal operating point: " + error.message() +
                " Loss loop at iterate " + std::to_string(iteration) + " (core losses at " + format_significant(coreTemperature) +
                " C, winding losses at " + format_significant(windingTemperature) + " C" +
                (coreHeldAtCeiling ? ", core held at its data limit" : "") + "). History:" + history,
                lossesAtAmbient.value());
        }
        catch (const CalculationException& error) {
            // Not handled: the same error, with the loss loop's iterates up to it, so the losses and
            // temperatures that drove the network there can be read.
            if (!lossTemperatureLoop) {
                throw CalculationException(error.code(), error.message() + " One-pass temperature estimate (loss-temperature "
                    "loop off): core and winding losses at the operating temperature, " + std::to_string(coreTemperature) + " C.");
            }
            throw CalculationException(error.code(), error.message() + " Loss loop at iterate " +
                std::to_string(iteration) + " (core losses at " + std::to_string(coreTemperature) + " C, winding losses at " +
                std::to_string(windingTemperature) + " C" + (coreHeldAtCeiling ? ", core held at its data limit" : "") +
                "). History:" + history);
        }
        if (!thermalResult.maximumCoreTemperature) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT,
                "Thermal network of '" + magnetic.get_reference() + "' has no core node: no core temperature to "
                "evaluate the core losses at");
        }
        if (withWinding && !thermalResult.averageTurnTemperature) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT,
                "Thermal network of '" + magnetic.get_reference() + "' has no turn node: no winding temperature to "
                "evaluate the winding losses at");
        }
        const double nextCoreTemperature = thermalResult.maximumCoreTemperature.value();
        const double nextWindingTemperature = withWinding ? thermalResult.averageTurnTemperature.value() : ambientTemperature;
        if (iteration == 1) {
            residualAtAmbient = std::make_pair(nextCoreTemperature - coreTemperature, nextWindingTemperature);
        }
        // Numbers a person reads: significant figures, scaled units (ABT #2001).
        history += " [" + std::to_string(iteration) + "] Pcore=" + format_quantity(coreLossesOutput.get_core_losses(), "W") +
                   " at " + format_significant(coreTemperature) + " C";
        if (magneticFluxDensity.get_processed() && magneticFluxDensity.get_processed()->get_peak()) {
            history += " (Bpk " + format_quantity(magneticFluxDensity.get_processed()->get_peak().value(), "T") + ")";
        }
        if (withWinding) {
            history += ", Pcu=" + format_quantity(windingLosses->get_winding_losses(), "W") + " at " +
                       format_significant(windingTemperature) + " C";
        }
        history += " -> core " + format_significant(nextCoreTemperature) + " C";
        if (withWinding) {
            history += ", winding " + format_significant(nextWindingTemperature) + " C";
        }
        history += ";";

        // One pass (loss-temperature loop off): these losses, at the operating temperature, and the hottest point
        // the network gives for them. Nothing is iterated, held or checked for runaway: there is no loop.
        // Held at the ceiling, the core's temperature is no longer iterated: only the winding's must settle.
        if (!lossTemperatureLoop ||
            ((coreHeldAtCeiling || std::fabs(nextCoreTemperature - coreTemperature) < kOperatingTemperatureTolerance) &&
             std::fabs(nextWindingTemperature - windingTemperature) < kOperatingTemperatureTolerance)) {
            // Converged (loop on): the losses were evaluated at temperatures within the tolerance of the
            // ones these very losses produce.
            SolvedLosses solved;
            solved.coreLosses = coreLossesOutput;
            solved.excitation = excitation;
            solved.windingLosses = windingLosses;
            solved.dcBiasPoint = dcBiasPoint;
            solved.lowestFixedPointSearchEvaluations = lowestFixedPointSearchEvaluations;
            solved.temperature.set_initial_temperature(config.ambientTemperature);
            // ABT #2078 (Alf, 2026-10-25): a core held at its data limit is at least as hot as the temperature it is
            // held at, whatever the network gives for the losses evaluated there. Held paramagnetic at Curie it loses
            // almost nothing and the network came back cold (a P 18/11 TM63 resonant inductor held at 289 C published
            // a 33 C hottest point), so a consumer reading only the temperature saw a cold, fine core. The published
            // hottest point is therefore never below the material's data limit the core is held at, and the method says
            // it is a floor.
            // The floor is the material's own limit (its Curie point, or where its data ends), never a number of ours.
            if (coreHeldAtCeiling && !dataCeiling) {
                throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT,
                                           "Loss-temperature solve: the core is held but its material has no data limit");
            }
            const double maximumTemperature = coreHeldAtCeiling ? std::max(thermalResult.maximumTemperature, dataCeiling->temperature)
                                                                : thermalResult.maximumTemperature;
            solved.temperature.set_maximum_temperature(maximumTemperature);
            double totalLosses = config.coreLosses + config.windingLosses;
            if (totalLosses > 0) {
                solved.temperature.set_bulk_thermal_resistance(
                    (thermalResult.maximumTemperature - config.ambientTemperature) / totalLosses);
            }
            solved.temperature.set_method_used(coreHeldAtCeiling ? "ThermalNetwork, core held at its material's data limit: "
                                                                   "the real temperature is at least this, and higher"
                                                                 : "ThermalNetwork");
            solved.temperature.set_origin(ResultOrigin::SIMULATION);
            return solved;
        }
        // Winding runaway (Alf, 2026-10-14, ABT #1855). Held at the ceiling, the core's losses and
        // temperature are fixed, and what is left is the one-dimensional map T_w -> g(T_w): the average
        // turn temperature the network returns for the winding losses evaluated at T_w. Its iteration
        // reaches a fixed point only where the map contracts, |g'| < 1. Two iterates of that same map give
        // its secant slope; when the winding is pushed hotter (g(T_w) > T_w), its losses grew between the
        // two iterates and the slope is at or above 1, the residual g(T_w) - T_w does not shrink as T_w
        // rises: the copper's losses heat it faster than the network sheds the heat, the iteration has no
        // fixed point ahead to settle on and the winding runs away (the 220 kW transformer: slope ~1.8,
        // 2.8 kW at 60 C, 2.3 MW and a 549,000 C hot spot fourteen iterates later). That is refused here,
        // from the iterates themselves, instead of iterating on until a limit or the network gives up.
        if (withWinding && coreHeldAtCeiling) {
            const double windingLossesAtThisIterate = windingLosses->get_winding_losses();
            if (heldWindingIterate) {
                const double step = windingTemperature - heldWindingIterate->windingTemperature;
                if (step > 0 && nextWindingTemperature > windingTemperature &&
                    windingLossesAtThisIterate > heldWindingIterate->windingLosses) {
                    const double slope = (nextWindingTemperature - heldWindingIterate->nextWindingTemperature) / step;
                    if (slope >= 1) {
                        std::ostringstream reason;
                        reason << "No thermal operating point: the winding runs away. '" << magnetic.get_reference()
                               << "' has its core held at " << format_significant(coreTemperature) << " C, the last temperature its material's data "
                               << "covers, and its winding losses grow from " << format_quantity(heldWindingIterate->windingLosses, "W") << " at "
                               << format_significant(heldWindingIterate->windingTemperature) << " C to " << format_quantity(windingLossesAtThisIterate, "W") << " at "
                               << format_significant(windingTemperature) << " C while the average turn temperature they produce grows from "
                               << format_significant(heldWindingIterate->nextWindingTemperature) << " C to " << format_significant(nextWindingTemperature)
                               << " C: dT_turn/dT_w = " << format_significant(slope) << ", at or above 1, where a fixed point needs it below 1."
                               << " History:" << history;
                        throw NoThermalOperatingPointException(ErrorCode::CALCULATION_DIVERGED, reason.str(),
                                                               lossesAtAmbient.value());
                    }
                }
            }
            heldWindingIterate = HeldWindingIterate{windingTemperature, nextWindingTemperature, windingLossesAtThisIterate};
        }

        // Next iterate. Plain substitution (relaxation 1) is the fastest when the losses' feedback
        // on temperature is weak, but a ferrite whose losses FALL steeply with temperature makes
        // the map overshoot the fixed point and swing around it forever (seen: 3.9 W at 106 C,
        // 2.6 W at 134 C, 0.01 W at 150 C). Each time a temperature's residual changes sign, that
        // temperature's step is halved, which turns the swing into convergence. Core and winding
        // keep their own step: damping the winding because the CORE swung would leave it crawling
        // toward its own answer. The convergence test above is on the residuals themselves, so
        // damping never makes an unsettled state look settled.
        const double coreResidual = nextCoreTemperature - coreTemperature;
        const double windingResidual = nextWindingTemperature - windingTemperature;
        if (previousCoreResidual && previousCoreResidual.value() * coreResidual < 0) {
            coreRelaxation *= 0.5;
        }
        if (previousWindingResidual && previousWindingResidual.value() * windingResidual < 0) {
            windingRelaxation *= 0.5;
        }
        previousCoreResidual = coreResidual;
        previousWindingResidual = windingResidual;
        double proposedCoreTemperature = coreTemperature + coreRelaxation * coreResidual;
        double proposedWindingTemperature = windingTemperature + windingRelaxation * windingResidual;

        // The core-loss and permeability models describe the FERROMAGNETIC material and refuse a
        // temperature at or above its Curie point; the permeability also refuses one between the end
        // of its measured mu_i(T) table and the Curie point, the magnetisation curve one where H_sat(T)
        // is not positive, and the Steinmetz losses one where ct(T) is not (ABT #1908). The first steps start from
        // ambient-temperature losses and can overshoot far past where the answer lies, so the iterate
        // is kept inside the range where the losses and the permeability exist: at most halfway from
        // the current core temperature to that ceiling.
        //
        // If the iteration has been driven to within the tolerance of the ceiling and the losses there
        // still push the core hotter, no operating point exists below it. Until ABT #1855 that was
        // refused. EXPLICIT, USER-APPROVED EXCEPTION to the no-fallback rule (Alf, 2026-10-14, ABT #1855):
        // the core is then held at the last temperature its material's data covers (this iterate, within
        // the tolerance of the ceiling), its core losses are evaluated there, and the winding temperature
        // is iterated on to its own fixed point with those core losses. Its permeability (Alf, 2026-10-14):
        //   - ceiling = Curie temperature: the core is magnetized as the paramagnetic core it is above
        //     its Curie point, mu_r = 1 (make_paramagnetic_magnetic);
        //   - ceiling = end of a mu_i(T) table that stops below Curie (3C96: 180 C, Curie above 200 C):
        //     the material is still magnetic there, so the core is held at the table's last temperature
        //     exactly and keeps the last tabulated permeability;
        //   - ceiling = where an evaluated Steinmetz range's ct(T) stops being positive (ABT #1908, the same hold as the mu_i(T) table's end): the
        //     material is still magnetic, and the models refuse at that temperature itself, so the core
        //     stays at this iterate, within the tolerance below it, with its own permeability there.
        // The numbers are therefore those of a core at its data limit: its real temperature is higher,
        // and ResultsAnalyzer::analyze_temperature says so for this operating point
        // (is_core_held_at_data_ceiling).
        if (coreHeldAtCeiling) {
            windingTemperature = proposedWindingTemperature;
            continue;
        }
        if (dataCeiling && proposedCoreTemperature >= dataCeiling->temperature) {
            const double distanceToCeiling = dataCeiling->temperature - coreTemperature;
            if (distanceToCeiling < kOperatingTemperatureTolerance) {
                // ABT #1840 (Alf, 2026-10-18: "find lowest crossing"). A part heated from ambient stops at the
                // LOWEST temperature where its losses and the network balance stably, and the iteration can miss
                // it: its first step from the ambient-temperature losses can land past an unstable fixed point
                // and climb to here (ETD 29 N87 at 0.44 of B_sat: g(25 C) = 186 C, stable point 97 C, unstable
                // 172 C). Before holding the core at its data limit, the bracket [ambient, ceiling] is searched
                // once for that point (find_lowest_stable_core_temperature below); the iteration resumes from it
                // and its own convergence test decides. With no fixed point the core is held, as before.
                if (!lowestFixedPointSearched) {
                    lowestFixedPointSearched = true;
                    auto lowestStable = findLowestStableCoreTemperature(ambientTemperature, coreTemperature,
                                                                        residualAtAmbient.value(), iteration);
                    if (lowestStable) {
                        history += " [search] lowest stable fixed point near " + format_significant(lowestStable->first) +
                                   " C after " + std::to_string(lowestFixedPointSearchEvaluations) + " evaluations, iteration resumes there;";
                        coreTemperature = lowestStable->first;
                        windingTemperature = lowestStable->second;
                        coreRelaxation = 1.0;
                        windingRelaxation = 1.0;
                        previousCoreResidual.reset();
                        previousWindingResidual.reset();
                        continue;
                    }
                    history += " [search] no fixed point below the data limit after " +
                               std::to_string(lowestFixedPointSearchEvaluations) + " evaluations;";
                }
                coreHeldAtCeiling = true;
                switch (dataCeiling->source) {
                    case CoreDataCeilingSource::CURIE_TEMPERATURE:
                        paramagneticMagnetic = make_paramagnetic_magnetic(magnetic);
                        break;
                    case CoreDataCeilingSource::PERMEABILITY_TABLE_END:
                        coreTemperature = dataCeiling->temperature;
                        break;
                    case CoreDataCeilingSource::STEINMETZ_TEMPERATURE_FACTOR_END:
                        // ABT #1908: at that temperature itself ct is zero and the model refuses,
                        // so the core stays at this iterate, within the tolerance below it, the material
                        // still magnetic and with its own permeability there.
                        break;
                }
                windingTemperature = proposedWindingTemperature;
                continue;
            }
            proposedCoreTemperature = coreTemperature + 0.5 * distanceToCeiling;
        }
        coreTemperature = proposedCoreTemperature;
        windingTemperature = proposedWindingTemperature;
    }
    // ABT #1960: the loop found no temperatures within its tolerance of the ones their losses produce. That is
    // the same answer as the network failing on an iterate's losses or the held winding running away: this
    // point has no thermal operating point the loop can report, and it is that point's finding, under the
    // caller's NoThermalOperatingPointPolicy (simulate: losses at ambient, no temperature; the advisers: THROW),
    // not a refusal of the whole design. Seen: a saturating 78 EL 25/6.3 inductor whose B(T) drops 0.637 ->
    // 0.630 T between 135.890 and 135.895 C, so its core losses jump 4.19 -> 3.81 W and the core temperature
    // they produce 138.1 -> 131.1 C: the iterates close in on that step from both sides and never land within
    // the tolerance.
    throw NoThermalOperatingPointException(ErrorCode::CALCULATION_DIVERGED,
        "No thermal operating point: the losses of '" + magnetic.get_reference() + "' and the temperatures they "
        "produce did not settle within " + std::to_string(kOperatingTemperatureMaximumIterations) +
        " iterations (tolerance " + format_significant(kOperatingTemperatureTolerance) + " K): no temperature was found "
        "within the tolerance of the one its own losses produce (thermal runaway, or losses that jump with temperature "
        "so the iterates close in on the jump from both sides). History:" + history,
        lossesAtAmbient.value());
}

namespace {
    // Builds a single-valued DimensionWithTolerance (nominal only) tagged with a unit.
    DimensionWithTolerance nominal_dimension(double value, const std::string& unit) {
        DimensionWithTolerance dimension;
        dimension.set_nominal(value);
        dimension.set_unit(unit);
        return dimension;
    }

    // RMS of a winding's current/voltage signal in an operating point, if present.
    std::optional<double> excitation_rms(const SignalDescriptor& signal) {
        if (!signal.get_processed() || !signal.get_processed()->get_rms()) {
            return std::nullopt;
        }
        return signal.get_processed()->get_rms();
    }
}

MagneticManufacturerInfo MagneticSimulator::build_datasheet(Mas& mas) {
    auto& magnetic = mas.get_mutable_magnetic();
    const auto& outputs = mas.get_outputs();
    const auto& operatingPoints = mas.get_inputs().get_operating_points();

    if (outputs.empty()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA,
            "Cannot build datasheet: MAS carries no simulation outputs; run MagneticSimulator::simulate first");
    }
    if (operatingPoints.size() != outputs.size()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA,
            "Cannot build datasheet: number of operating points (" + std::to_string(operatingPoints.size()) +
            ") does not match number of outputs (" + std::to_string(outputs.size()) + ")");
    }

    const auto& windings = magnetic.get_coil().get_functional_description();
    size_t numberWindings = windings.size();

    // ----------------------------------------------------------------- Part
    Part part;
    part.set_number_of_windings(static_cast<int64_t>(numberWindings));
    auto materialName = magnetic.get_mutable_core().get_material_name();
    if (!materialName.empty()) {
        part.set_material(materialName);
    }
    if (magnetic.get_manufacturer_info()) {
        auto existing = magnetic.get_manufacturer_info().value();
        if (existing.get_reference())   part.set_part_number(existing.get_reference());
        if (existing.get_description()) part.set_description(existing.get_description());
    }

    // ----------------------------------------------------------- Electrical
    MagneticDatasheetElectrical electrical;
    electrical.set_subtype(numberWindings > 1 ? ElectricalSubtype::COUPLED_INDUCTOR : ElectricalSubtype::INDUCTOR);

    // Magnetizing inductance and leakage inductance are taken from the nominal (first) operating point.
    if (outputs[0].get_inductance()) {
        auto inductanceOutput = outputs[0].get_inductance().value();
        electrical.set_inductance(inductanceOutput.get_magnetizing_inductance().get_magnetizing_inductance());

        // Per-winding leakage array is winding-indexed with 0 at the primary slot, so the
        // first secondary lives at index 1.
        if (inductanceOutput.get_leakage_inductance() &&
            inductanceOutput.get_leakage_inductance()->get_leakage_inductance_per_winding().size() > 1) {
            electrical.set_leakage_inductance(
                inductanceOutput.get_leakage_inductance()->get_leakage_inductance_per_winding()[1]);
        }
    }

    // DC resistance per winding, evaluated at the nominal operating point's ambient temperature, via MKF.
    double dcResistanceTemperature = operatingPoints[0].get_conditions().get_ambient_temperature();
    auto dcResistancePerWinding = WindingOhmicLosses::calculate_dc_resistance_per_winding(
        magnetic.get_coil(), dcResistanceTemperature);
    if (numberWindings == 1 && dcResistancePerWinding.size() == 1) {
        electrical.set_dc_resistance(nominal_dimension(dcResistancePerWinding[0], "Ohm"));
    }
    else if (dcResistancePerWinding.size() == numberWindings) {
        std::vector<DimensionWithTolerance> dcResistances;
        for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
            dcResistances.push_back(nominal_dimension(dcResistancePerWinding[windingIndex], "Ohm"));
        }
        electrical.set_dc_resistances(dcResistances);
    }

    // Rated current: the heat-limited DC current for the standard 40 K rise (datasheet Irms),
    // computed from the winding ohmic loss and the core temperature model.
    electrical.set_rated_currents(std::vector<double>{magnetic.calculate_rated_current()});

    // Self-resonant frequency and an impedance-vs-frequency sweep. The sweep uses the wideband
    // Foster-ladder model (the magnetizing resonance plus one leakage resonance per extra
    // winding), built once and evaluated per frequency; stray capacitance and winding
    // resistance are folded in by the model. Anchored on the SRF and spanning several decades
    // above it so the published curve captures the higher (leakage) resonances too.
    {
        Impedance impedanceModel;
        double selfResonantFrequency = impedanceModel.calculate_self_resonant_frequency(magnetic, dcResistanceTemperature);
        electrical.set_self_resonant_frequency(selfResonantFrequency);

        auto widebandModel = impedanceModel.build_wideband_impedance_model(magnetic, selfResonantFrequency, dcResistanceTemperature);

        const size_t numberPoints = 51;
        double startFrequency = selfResonantFrequency / 1000;
        double stopFrequency = selfResonantFrequency * 100;
        // The published curve covers only the frequencies where the core material has complex
        // permeability data: outside its table there is nothing to evaluate (get_complex_permeability
        // throws there). This sweep is the datasheet's own choice of band, so it asks for the
        // material's range explicitly and keeps the part of the band that lies inside it. If none
        // of the band does, the curve cannot be computed and that is reported, not filled in.
        {
            auto coreMaterial = magnetic.get_core().resolve_material();
            auto [minimumMaterialFrequency, maximumMaterialFrequency] = ComplexPermeability().get_frequency_range(coreMaterial);
            double requestedStartFrequency = startFrequency;
            double requestedStopFrequency = stopFrequency;
            startFrequency = std::max(startFrequency, minimumMaterialFrequency);
            stopFrequency = std::min(stopFrequency, maximumMaterialFrequency);
            if (!(startFrequency < stopFrequency)) {
                throw ComplexPermeabilityFrequencyOutOfRangeException(coreMaterial.get_name(),
                    requestedStartFrequency < minimumMaterialFrequency ? requestedStopFrequency : requestedStartFrequency,
                    minimumMaterialFrequency, maximumMaterialFrequency);
            }
        }
        // Exact ends (logarithmic_spaced_array): start * 10^(n * step) can land one ulp past the
        // material's last tabulated frequency, where complex permeability throws.
        auto frequencies = logarithmic_spaced_array(startFrequency, stopFrequency, numberPoints);

        std::vector<DatasheetImpedancePoint> impedancePoints;
        double maximumImpedanceMagnitude = std::numeric_limits<double>::lowest();
        for (size_t pointIndex = 0; pointIndex < numberPoints; ++pointIndex) {
            double frequency = frequencies[pointIndex];
            auto impedance = impedanceModel.impedance_from_model(widebandModel, frequency);

            ImpedancePoint impedancePoint;
            impedancePoint.set_magnitude(std::abs(impedance));
            impedancePoint.set_real_part(impedance.real());
            impedancePoint.set_imaginary_part(impedance.imag());
            impedancePoint.set_phase(std::arg(impedance));

            DatasheetImpedancePoint datasheetImpedancePoint;
            datasheetImpedancePoint.set_frequency(frequency);
            datasheetImpedancePoint.set_impedance(impedancePoint);
            impedancePoints.push_back(datasheetImpedancePoint);

            maximumImpedanceMagnitude = std::max(maximumImpedanceMagnitude, std::abs(impedance));
        }
        electrical.set_impedance_points(impedancePoints);
        electrical.set_maximum_impedance(maximumImpedanceMagnitude);
    }

    // Turns ratio (primary : first secondary) for coupled parts.
    if (numberWindings > 1) {
        auto turnsRatios = magnetic.get_turns_ratios();
        if (!turnsRatios.empty()) {
            electrical.set_turns_ratios(std::vector<DimensionWithTolerance>{nominal_dimension(turnsRatios[0], "")});
        }
    }

    // Peak saturation current, evaluated at the hottest operating point (lowest, conservative I_sat),
    // delegating to the MKF saturation-current model.
    size_t hottestOperatingPointIndex = 0;
    double hottestAmbient = std::numeric_limits<double>::lowest();
    for (size_t operatingPointIndex = 0; operatingPointIndex < operatingPoints.size(); ++operatingPointIndex) {
        double ambient = operatingPoints[operatingPointIndex].get_conditions().get_ambient_temperature();
        if (ambient > hottestAmbient) {
            hottestAmbient = ambient;
            hottestOperatingPointIndex = operatingPointIndex;
        }
    }
    {
        OperatingPoint saturationOperatingPoint = operatingPoints[hottestOperatingPointIndex];
        electrical.set_saturation_current_peak(
            magnetic.calculate_saturation_current(saturationOperatingPoint, hottestAmbient));
    }

    // -------------------------------------------------------------- Thermal
    Thermal thermal;

    // Operating temperature range = span of ambient conditions across operating points.
    double minimumAmbient = std::numeric_limits<double>::max();
    double maximumAmbient = std::numeric_limits<double>::lowest();
    for (const auto& operatingPoint : operatingPoints) {
        double ambient = operatingPoint.get_conditions().get_ambient_temperature();
        minimumAmbient = std::min(minimumAmbient, ambient);
        maximumAmbient = std::max(maximumAmbient, ambient);
    }
    {
        DimensionWithTolerance operatingTemperature;
        operatingTemperature.set_minimum(minimumAmbient);
        operatingTemperature.set_maximum(maximumAmbient);
        operatingTemperature.set_nominal(operatingPoints[0].get_conditions().get_ambient_temperature());
        operatingTemperature.set_unit("Celsius");
        thermal.set_operating_temperature(operatingTemperature);
    }

    // Worst-case hot-spot temperature rise and thermal resistance, computed with the full
    // Temperature thermal-network model from the simulated core + winding losses at each
    // operating point (the per-turn winding-loss distribution is passed through for fidelity).
    Magnetic thermalMagnetic = magnetic;
    if (!thermalMagnetic.get_coil().get_turns_description()) {
        thermalMagnetic.get_mutable_coil().wind();
    }
    double worstTemperatureRise = std::numeric_limits<double>::lowest();
    std::optional<double> worstThermalResistance;
    for (size_t operatingPointIndex = 0; operatingPointIndex < outputs.size(); ++operatingPointIndex) {
        if (!outputs[operatingPointIndex].get_core_losses()) {
            continue;
        }
        double ambient = operatingPoints[operatingPointIndex].get_conditions().get_ambient_temperature();
        // Shared config builder (ABT #906) — same as outputs[].temperature and the plot wrappers.
        auto config = TemperatureConfig::fromSimulatedOutput(operatingPoints[operatingPointIndex],
                                                             outputs[operatingPointIndex]);

        double hotspot = Temperature(thermalMagnetic, config).calculateTemperatures().maximumTemperature;
        double temperatureRise = hotspot - ambient;
        if (temperatureRise > worstTemperatureRise) {
            worstTemperatureRise = temperatureRise;
            double totalLosses = config.coreLosses + config.windingLosses;
            worstThermalResistance = totalLosses > 0 ? std::optional<double>(temperatureRise / totalLosses)
                                                     : std::nullopt;
        }
    }
    if (worstTemperatureRise > std::numeric_limits<double>::lowest()) {
        thermal.set_temperature_rise(worstTemperatureRise);
        if (worstThermalResistance) {
            thermal.set_thermal_resistance(worstThermalResistance.value());
        }
    }

    // ----------------------------------------------------------- Mechanical
    Mechanical mechanical;
    auto dimensions = magnetic.get_maximum_dimensions();  // [width, height, depth] in metres
    mechanical.set_width(nominal_dimension(dimensions[0], "m"));
    mechanical.set_height(nominal_dimension(dimensions[1], "m"));
    mechanical.set_length(nominal_dimension(dimensions[2], "m"));

    // ---------------------------------------------------------- Application
    // The application block states the circuit conditions this design was simulated for.
    MagneticDatasheetApplication application;
    application.set_switching_frequency(operatingPoints[0].get_excitations_per_winding()[0].get_frequency());

    // Input voltage range (primary winding) across operating points.
    double minimumInputVoltage = std::numeric_limits<double>::max();
    double maximumInputVoltage = std::numeric_limits<double>::lowest();
    for (const auto& operatingPoint : operatingPoints) {
        const auto& excitation = operatingPoint.get_excitations_per_winding()[0];
        if (excitation.get_voltage()) {
            auto rms = excitation_rms(excitation.get_voltage().value());
            if (rms) {
                minimumInputVoltage = std::min(minimumInputVoltage, rms.value());
                maximumInputVoltage = std::max(maximumInputVoltage, rms.value());
            }
        }
    }
    if (maximumInputVoltage > std::numeric_limits<double>::lowest()) {
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_minimum(minimumInputVoltage);
        inputVoltage.set_maximum(maximumInputVoltage);
        inputVoltage.set_unit("V");
        application.set_input_voltage(Dimension(inputVoltage));
    }

    // Output voltages (secondary windings) at the nominal operating point.
    if (numberWindings > 1 &&
        operatingPoints[0].get_excitations_per_winding().size() == numberWindings) {
        std::vector<double> outputVoltages;
        for (size_t windingIndex = 1; windingIndex < numberWindings; ++windingIndex) {
            const auto& excitation = operatingPoints[0].get_excitations_per_winding()[windingIndex];
            if (excitation.get_voltage()) {
                auto rms = excitation_rms(excitation.get_voltage().value());
                if (rms) {
                    outputVoltages.push_back(rms.value());
                }
            }
        }
        if (!outputVoltages.empty()) {
            application.set_output_voltages(outputVoltages);
        }
    }

    // ------------------------------------------------- Assemble datasheetInfo
    DatasheetInfo datasheetInfo;
    datasheetInfo.set_part(part);
    datasheetInfo.set_electrical(std::vector<MagneticDatasheetElectrical>{electrical});
    datasheetInfo.set_thermal(thermal);
    datasheetInfo.set_mechanical(mechanical);
    datasheetInfo.set_application(application);

    // Preserve any existing manufacturer metadata; only (re)populate the datasheetInfo block.
    MagneticManufacturerInfo manufacturerInfo;
    if (magnetic.get_manufacturer_info()) {
        manufacturerInfo = magnetic.get_manufacturer_info().value();
    }
    else {
        manufacturerInfo.set_name("OpenMagnetics");
    }
    manufacturerInfo.set_datasheet_info(datasheetInfo);

    magnetic.set_manufacturer_info(manufacturerInfo);
    return manufacturerInfo;
}

} // namespace OpenMagnetics
