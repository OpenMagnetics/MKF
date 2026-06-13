#include "processors/Inputs.h"
#include <source_location>
#include <cfloat>
#include <cmath>
#include <limits>
#include "advisers/MagneticAdviser.h"
#include "advisers/CoreAdviser.h"
#include "advisers/MagneticFilterInternal.h"  // is_energy_storing_topology()
#include "physical_models/Impedance.h"
#include "processors/MagneticSimulator.h"
#include "support/Painter.h"
#include <magic_enum_utility.hpp>
#include "support/Logger.h"


namespace OpenMagnetics {

namespace {
// RAII: back up the global core/wire databases, replace them with constraint-
// filtered copies for the duration of the scope, and restore on destruction
// (also on exception). A no-op when the relevant section of `constraints` is
// empty — safe to construct unconditionally at the top of any ctx-aware
// adviser overload.
struct DatabaseFilterScope {
    explicit DatabaseFilterScope(const AdviserConstraints& constraints) {
        if (!constraints.shapeFamily.empty() || !constraints.coreMaterialType.empty()) {
            if (coreDatabase.empty()) load_cores();
            _coreBackup = coreDatabase;
            coreDatabase = filterCoresByConstraints(coreDatabase, constraints);
            _restoreCores = true;
        }
        if (!constraints.wireType.empty()) {
            if (wireDatabase.empty()) load_wires();
            _wireBackup = wireDatabase;
            wireDatabase = filterWiresByConstraints(wireDatabase, constraints);
            _restoreWires = true;
        }
    }
    ~DatabaseFilterScope() {
        if (_restoreCores) coreDatabase = std::move(_coreBackup);
        if (_restoreWires) wireDatabase = std::move(_wireBackup);
    }
    DatabaseFilterScope(const DatabaseFilterScope&) = delete;
    DatabaseFilterScope& operator=(const DatabaseFilterScope&) = delete;
private:
    std::vector<Core> _coreBackup;
    std::map<std::string, Wire> _wireBackup;
    bool _restoreCores = false;
    bool _restoreWires = false;
};
} // namespace

void MagneticAdviser::set_unique_core_shapes(bool value) {
    _uniqueCoreShapes = value;
}

bool MagneticAdviser::get_unique_core_shapes() {
    return _uniqueCoreShapes;
}

void MagneticAdviser::set_application(Application value) {
    _application = value;
}

Application MagneticAdviser::get_application() {
    return _application;
}

void MagneticAdviser::set_core_mode(CoreAdviser::CoreAdviserModes value) {
    _coreAdviserMode = value;
}

CoreAdviser::CoreAdviserModes MagneticAdviser::get_core_mode() {
    return _coreAdviserMode;
}

void MagneticAdviser::load_filter_flow(std::vector<MagneticFilterOperation> flow, std::optional<Inputs> inputs) {
    _filters.clear();
    _loadedFilterFlow = flow;
    for (auto filterConfiguration : flow) {
        MagneticFilters filterEnum = filterConfiguration.get_filter();
        _filters[filterEnum] = MagneticFilter::factory(filterEnum, inputs);
    }
}

std::vector<std::pair<Mas, double>> MagneticAdviser::get_advised_magnetic_fast(Inputs inputs, size_t maximumNumberResults) {
    inputs = pre_process_inputs(inputs);

    if (coreDatabase.empty()) {
        load_cores();
    }

    CoreAdviser coreAdviser;
    coreAdviser.set_application(get_application());
    coreAdviser.set_mode(get_core_mode());

    // Step 1: Create Magnetic objects from all cores in the database
    auto magneticsWithScoring = coreAdviser.create_magnetic_dataset(inputs, &coreDatabase, false);

    // Step 2: Area product filter — sorted by AP score, fast binary filter
    CoreAdviser::MagneticCoreFilterAreaProduct filterAreaProduct(inputs);
    filterAreaProduct.set_scorings(&coreAdviser._scorings);
    filterAreaProduct.set_filter_configuration(&coreAdviser._filterConfiguration);
    magneticsWithScoring = filterAreaProduct.filter_magnetics(&magneticsWithScoring, inputs, 1.0, true);

    if (magneticsWithScoring.empty()) {
        return {};
    }

    // Step 2b: Energy stored filter — reject cores that cannot store the required energy.
    // Without this, tiny cores pass the AP check but produce absurd losses for
    // energy-storing topologies (e.g. flyback) because they need impossibly large gaps.
    std::map<std::string, std::string> defaultModels{
        {"gapReluctance", to_string(defaults.reluctanceModelDefault)},
        {"coreLosses", to_string(defaults.coreLossesModelDefault)},
        {"coreTemperature", to_string(defaults.coreTemperatureModelDefault)}
    };
    CoreAdviser::MagneticCoreFilterEnergyStored filterEnergyStored(inputs, defaultModels);
    filterEnergyStored.set_scorings(&coreAdviser._scorings);
    filterEnergyStored.set_filter_configuration(&coreAdviser._filterConfiguration);
    magneticsWithScoring = filterEnergyStored.filter_magnetics(&magneticsWithScoring, inputs, 1.0, true);

    if (magneticsWithScoring.empty()) {
        return {};
    }

    // Cap candidates: AP filter sorts by score, keep the best ones.
    // post_process_core takes ~2-3ms each; 100 candidates ≈ 250ms.
    const size_t maxCandidates = std::max(maximumNumberResults * 20, size_t(50));
    if (magneticsWithScoring.size() > maxCandidates) {
        magneticsWithScoring.resize(maxCandidates);
    }

    // Step 2c: For energy-storing topologies (flyback, buck, boost,
    // SEPIC, Cuk, Zeta, isolated-buck, PFC, CMC/DMC, isolated buck-
    // boost), toroidal ferrite cores are physically incompatible
    // with the core's job: ferrite has no distributed gap (unlike
    // powder), and a toroid has no parting plane for a discrete air
    // gap. Result: ungapped high-µ ferrite saturates at a fraction
    // of an ampere — useless for flyback-class energy storage.
    //
    // Without this filter, the fast path served ferrite toroids for
    // flyback designs (Heaviside corpus: T 38.1/19.05/12.7 with N30,
    // isat=1.2 A vs ipeak=1.8 A — see docs/MKF-HANDOFF.md). The
    // slow-path CoreAdviserPipeline already splits toroidal->powder
    // and non-toroidal->ferrite paths so this combination never
    // arises; the fast path needed the same guard.
    //
    // We do NOT call add_gapping_standard_cores here (the slow path
    // does). Non-toroidal ferrite cores are left ungapped at this
    // stage; add_initial_turns_by_inductance computes turns against
    // whatever gap the core has, and filterSaturation rejects any
    // candidate whose B_peak exceeds B_sat. That's the existing
    // fast-path contract — this step only narrows it physically.
    {
        auto topology = inputs.get_design_requirements().get_topology();
        if (is_energy_storing_topology(topology)) {
            std::vector<std::pair<Magnetic, double>> filtered;
            filtered.reserve(magneticsWithScoring.size());
            for (auto& entry : magneticsWithScoring) {
                const auto& core = entry.first.get_core();
                bool isToroidal = core.get_functional_description().get_type()
                                  == CoreType::TOROIDAL;
                if (isToroidal) {
                    auto material = core.resolve_material();
                    if (material.get_material() == MaterialType::FERRITE) {
                        // Drop — ungappable ferrite toroid for energy-storing app.
                        continue;
                    }
                }
                filtered.push_back(std::move(entry));
            }
            magneticsWithScoring = std::move(filtered);

            if (magneticsWithScoring.empty()) {
                return {};
            }
        }
    }

    // Step 2d: Gap non-toroidal ferrite cores. Without a gap, high-µ
    // ferrite achieves the target L with very few turns and B_peak
    // immediately saturates the core — useless for any application
    // that stores meaningful flux. Matches the slow path
    // (CoreAdviserPipeline::add_gapping_standard_cores call at
    // CoreAdviserPipeline.cpp:578, before add_initial_turns_by_inductance).
    // Toroidal cores (powder) are skipped: their distributed gap is
    // already baked into μ_eff.
    {
        std::vector<std::pair<Magnetic, double>> nonToroidal;
        std::vector<std::pair<Magnetic, double>> toroidal;
        nonToroidal.reserve(magneticsWithScoring.size());
        for (auto& entry : magneticsWithScoring) {
            if (entry.first.get_core().get_functional_description().get_type()
                == CoreType::TOROIDAL) {
                toroidal.push_back(std::move(entry));
            } else {
                nonToroidal.push_back(std::move(entry));
            }
        }
        if (!nonToroidal.empty()) {
            coreAdviser.add_gapping_standard_cores(&nonToroidal, inputs);
        }
        magneticsWithScoring = std::move(nonToroidal);
        for (auto& entry : toroidal) {
            magneticsWithScoring.push_back(std::move(entry));
        }
    }

    // Step 3: Set turns and gap analytically (single pass, no iteration)
    add_initial_turns_by_inductance(&magneticsWithScoring, inputs);

    // Step 3a: Inductance-validity filter — reject cores whose ACHIEVED
    // inductance is outside the required tolerance band. This is the slow
    // path's MagneticCoreFilterMagneticInductance, which the fast path was
    // missing. Without it, integer-turns rounding on small high-frequency
    // inductors can overshoot the target L several-fold (e.g. 2.4µH target
    // → 14µH actual); the over-inductance core then PASSES the saturation
    // filter precisely because its (wrong) high L gives low ripple and low
    // peak current. Saturation cannot catch wrong-inductance — only this
    // filter can. Must run after turns are set, before saturation.
    CoreAdviser::MagneticCoreFilterMagneticInductance filterMagneticInductance;
    filterMagneticInductance.set_scorings(&coreAdviser._scorings);
    filterMagneticInductance.set_filter_configuration(&coreAdviser._filterConfiguration);
    magneticsWithScoring = filterMagneticInductance.filter_magnetics(&magneticsWithScoring, inputs, 1, true);

    // Step 3b: Saturation filter — reject cores that exceed flux density saturation
    // after turns have been set. Must come after add_initial_turns_by_inductance.
    CoreAdviser::MagneticCoreFilterSaturation filterSaturation;
    filterSaturation.set_scorings(&coreAdviser._scorings);
    filterSaturation.set_filter_configuration(&coreAdviser._filterConfiguration);
    magneticsWithScoring = filterSaturation.filter_magnetics(&magneticsWithScoring, inputs, 1, true);

    // Step 4: Add secondary windings from turns ratios
    correct_windings(&magneticsWithScoring, inputs);

    // Step 5: Evaluate each core with fast_wind + ohmic losses + core losses
    std::vector<std::pair<Mas, double>> results;
    for (auto& [magnetic, scoring] : magneticsWithScoring) {
        try {
            auto mas = coreAdviser.post_process_core(magnetic, inputs);

            // Score by total losses (lower is better)
            double totalLosses = 0;
            for (auto& output : mas.get_outputs()) {
                if (output.get_core_losses()) {
                    totalLosses += output.get_core_losses()->get_core_losses();
                }
                if (output.get_winding_losses()) {
                    totalLosses += output.get_winding_losses()->get_winding_losses();
                }
            }

            if (totalLosses > 0 && std::isfinite(totalLosses)) {
                results.push_back({mas, totalLosses});
            }
        }
        catch (const std::exception& e) {
            logEntry(std::string("MagneticAdviser::get_advised_magnetic_fast: skipping candidate, scoring failed: ") + e.what(), "MagneticAdviser", 2);
            continue;
        }
    }

    // Step 6: Sort by total losses (ascending — lower losses = better)
    std::stable_sort(results.begin(), results.end(),
        [](const auto& a, const auto& b) { return a.second < b.second; });

    // Step 7: Truncate to requested number
    if (results.size() > maximumNumberResults) {
        results.resize(maximumNumberResults);
    }

    return results;
}

std::vector<std::pair<Mas, double>> MagneticAdviser::get_advised_magnetic(Inputs inputs, size_t maximumNumberResults) {
    return get_advised_magnetic(inputs, _defaultCustomMagneticFilterFlow, maximumNumberResults);
}

std::vector<std::pair<Mas, double>> MagneticAdviser::get_advised_magnetic(
        Inputs inputs,
        size_t maximumNumberResults,
        const LibraryContext* ctx,
        const AdviserConstraints& constraints) {
    auto scope = ctx ? ctx->applyScoped() : LibraryContext::Scope{};
    DatabaseFilterScope dbScope(constraints);
    return get_advised_magnetic(inputs, _defaultCustomMagneticFilterFlow, maximumNumberResults);
}

std::vector<std::pair<Mas, double>> MagneticAdviser::get_advised_magnetic(
        Inputs inputs,
        std::map<MagneticFilters, double> weights,
        size_t maximumNumberResults,
        const LibraryContext* ctx,
        const AdviserConstraints& constraints) {
    auto scope = ctx ? ctx->applyScoped() : LibraryContext::Scope{};
    DatabaseFilterScope dbScope(constraints);
    return get_advised_magnetic(inputs, weights, maximumNumberResults);
}

std::vector<std::pair<Mas, double>> MagneticAdviser::get_advised_magnetic(
        Inputs inputs,
        std::vector<MagneticFilterOperation> filterFlow,
        size_t maximumNumberResults,
        const LibraryContext* ctx,
        const AdviserConstraints& constraints) {
    auto scope = ctx ? ctx->applyScoped() : LibraryContext::Scope{};
    DatabaseFilterScope dbScope(constraints);
    return get_advised_magnetic(inputs, filterFlow, maximumNumberResults);
}

std::vector<std::pair<Mas, double>> MagneticAdviser::get_advised_magnetic_fast(
        Inputs inputs,
        size_t maximumNumberResults,
        const LibraryContext* ctx,
        const AdviserConstraints& constraints) {
    auto scope = ctx ? ctx->applyScoped() : LibraryContext::Scope{};
    DatabaseFilterScope dbScope(constraints);
    return get_advised_magnetic_fast(inputs, maximumNumberResults);
}

std::vector<std::pair<Mas, double>> MagneticAdviser::get_advised_magnetic(Inputs inputs, std::map<MagneticFilters, double> weights, size_t maximumNumberResults) {
    std::vector<MagneticFilterOperation> customMagneticFilterFlow{
        MagneticFilterOperation(MagneticFilters::COST, true, true, weights[MagneticFilters::COST]),
        MagneticFilterOperation(MagneticFilters::LOSSES, true, true, weights[MagneticFilters::LOSSES]),
        MagneticFilterOperation(MagneticFilters::DIMENSIONS, true, true, weights[MagneticFilters::DIMENSIONS]),
    };
    return get_advised_magnetic(inputs, customMagneticFilterFlow, maximumNumberResults);
}

std::vector<std::pair<Mas, double>> MagneticAdviser::get_advised_magnetic(Inputs inputs, std::vector<MagneticFilterOperation> filterFlow, size_t maximumNumberResults) {
    clear_scoring();
    load_filter_flow(filterFlow, inputs);
    std::vector<Mas> masData;

    // Adopt the application tag from the inputs when the caller didn't set one
    // explicitly. CMC / DMC wizards tag their DesignRequirements with
    // INTERFERENCE_SUPPRESSION so the downstream filters (toroidal-only,
    // CM-appropriate ferrites, bifilar winding) can kick in — without this
    // auto-adoption the adviser defaults to POWER and treats the CMC like a
    // generic transformer. INTERFERENCE_SUPPRESSION also implies
    // AVAILABLE_CORES mode: toroidal EMI cores are catalog parts you pick off
    // the shelf, not blanks you grind a gap into.
    {
        auto reqAppOpt = inputs.get_design_requirements().get_application();
        if (reqAppOpt.has_value() && _application == Application::POWER) {
            set_application(reqAppOpt.value());
            if (reqAppOpt.value() == Application::INTERFERENCE_SUPPRESSION &&
                _coreAdviserMode == CoreAdviser::CoreAdviserModes::STANDARD_CORES) {
                set_core_mode(CoreAdviser::CoreAdviserModes::AVAILABLE_CORES);
            }
        }
    }

    // Store original toroid setting for potential retry
    bool toroidsOriginallyEnabled = settings.get_use_toroidal_cores();

    // Snapshot all three core-filter settings so we can restore them on exit.
    // Without this, INTERFERENCE_SUPPRESSION (CMC/DMC) calls leak their
    // toroidal-only / concentric-disabled / cores-in-stock overrides into the
    // process-wide `settings` singleton AND into the cached global
    // `coreDatabase` (which is filtered at load time and only re-loaded when
    // empty). A subsequent caller — e.g. a DAB or LLC adviser test — then
    // sees an already-populated DB containing only toroids, prunes all
    // ferrite shapes away, and returns zero magnetics. See test pollution
    // bug: TestTopologyDabAdviser / TestTopologyLlc REQUIRE(results>=1) fails
    // when run after TestTopologyCmc.
    bool savedUseToroidalCores = settings.get_use_toroidal_cores();
    bool savedUseConcentricCores = settings.get_use_concentric_cores();
    bool savedUseOnlyCoresInStock = settings.get_use_only_cores_in_stock();

    // RAII restore of the three core-filter settings + the cached core DB. The
    // body below mutates these globals (high-voltage / interference-suppression
    // paths) and runs many throwing operations (get_advised_core, coil adviser,
    // simulate, score_magnetics); a manual restore at the end leaks the mutated
    // flags into every later call if any of them throws. The destructor restores
    // the snapshot on every exit path and drops the cached DB if a filter
    // actually changed, so the next consumer rebuilds it from the restored set.
    struct CoreFilterSettingsRestorer {
        Settings& settings;
        bool toroidal, concentric, onlyInStock;
        ~CoreFilterSettingsRestorer() {
            bool changed = settings.get_use_toroidal_cores() != toroidal ||
                           settings.get_use_concentric_cores() != concentric ||
                           settings.get_use_only_cores_in_stock() != onlyInStock;
            settings.set_use_toroidal_cores(toroidal);
            settings.set_use_concentric_cores(concentric);
            settings.set_use_only_cores_in_stock(onlyInStock);
            if (changed) {
                clear_loaded_cores();
            }
        }
    } coreFilterSettingsRestorer{settings, savedUseToroidalCores, savedUseConcentricCores, savedUseOnlyCoresInStock};

    // Check if high voltage requirements make toroids impractical
    // High voltage (>600V) with strict insulation requirements rarely works with toroids
    double maxVoltage = 0;
    auto insulationOpt = inputs.get_design_requirements().get_insulation();
    if (insulationOpt) {
        auto mainSupplyVoltageOpt = insulationOpt.value().get_main_supply_voltage();
        if (mainSupplyVoltageOpt) {
            auto voltages = mainSupplyVoltageOpt.value();
            maxVoltage = std::max({voltages.get_nominal().value_or(0), 
                                   voltages.get_minimum().value_or(0), 
                                   voltages.get_maximum().value_or(0)});
        }
    }
    
    if (toroidsOriginallyEnabled && maxVoltage > 600) {
        logEntry("High voltage requirements (" + std::to_string(int(maxVoltage)) + "V) detected. Disabling toroidal cores for better results.", "MagneticAdviser", 0);
        settings.set_use_toroidal_cores(false);
        toroidsOriginallyEnabled = false; // Don't retry since we proactively disabled
        clear_loaded_cores();
    }

    if (get_application() == Application::INTERFERENCE_SUPPRESSION) {
        settings.set_use_toroidal_cores(true);
        settings.set_use_only_cores_in_stock(false);
        settings.set_use_concentric_cores(false);
        // Force the global core DB to be rebuilt with the new filter settings
        // (it's filtered at load time and otherwise only reloaded when empty).
        clear_loaded_cores();
    }

    if (coreDatabase.empty()) {
        load_cores();
    }
    if (wireDatabase.empty()) {
        load_wires();
    }

    bool previousCoilIncludeAdditionalCoordinates = settings.get_coil_include_additional_coordinates();
    // RAII: keep the bool for the in-loop flag check below, but restore the
    // setting on every exit path (the body has throwing call paths).
    SettingsGuard<bool> coilIncludeCoordinatesGuard(settings,
        &Settings::get_coil_include_additional_coordinates,
        &Settings::set_coil_include_additional_coordinates, false);

    std::map<CoreAdviser::CoreAdviserFilters, double> coreWeights;
    for (auto flowStep : filterFlow) {
        if (flowStep.get_filter() == MagneticFilters::COST) {
            coreWeights[CoreAdviser::CoreAdviserFilters::COST] = flowStep.get_weight();
        }
        if (flowStep.get_filter() == MagneticFilters::DIMENSIONS) {
            coreWeights[CoreAdviser::CoreAdviserFilters::DIMENSIONS] = flowStep.get_weight();
        }
        if (flowStep.get_filter() == MagneticFilters::LOSSES) {
            coreWeights[CoreAdviser::CoreAdviserFilters::EFFICIENCY] = flowStep.get_weight();
        }
    }

    CoreAdviser coreAdviser;

    coreAdviser.set_unique_core_shapes(true);
    coreAdviser.set_application(get_application());
    coreAdviser.set_mode(get_core_mode());
    CoilAdviser coilAdviser;
    MagneticSimulator magneticSimulator;
    size_t numberWindings = inputs.get_design_requirements().get_turns_ratios().size() + 1;
    size_t coresWound = 0;

    logEntry("Getting core", "MagneticAdviser", 2);
    size_t expectedWoundCores = std::min(maximumNumberResults, std::max(size_t(2), size_t(floor(double(maximumNumberResults) / numberWindings))));
    size_t requestedCores = expectedWoundCores;
    std::vector<std::string> evaluatedCores;
    size_t previouslyObtainedCores = SIZE_MAX;
    size_t whileIteration = 0;
    const size_t maxWhileIterations = 2;  // Limit exponential growth
    const size_t maxEvaluatedCores = 40;  // OPTIMIZATION: Balanced - 40 cores to capture variety including medium-sized
    // Memory bounds for the candidate pool. Each Mas in `masData` carries a
    // full ngspice-simulated operating point (per-turn time-series waveforms),
    // which can run to MBs per candidate. Without these caps, with
    // maximumNumberResults=100 we used to retain up to 40 cores * 50 coils
    // per core = 2000 fully-simulated mas, exhausting the WASM heap
    // (std::bad_alloc). See MKF_ISSUES.md M3.
    //   - Per-core cap: at most 5 coils per core, while still respecting very
    //     small request sizes (min with the original 0.5*N expression).
    //   - Global cap: never grow the pool past maximumNumberResults*4 — final
    //     scoring only keeps the top maximumNumberResults anyway, so 4x gives
    //     ample variety for the scoring/sort step without unbounded growth.
    const size_t perCoreCoilCap = std::min(size_t(5), size_t(ceil(maximumNumberResults * 0.5)));
    const size_t globalCandidateCap = std::max(size_t(1), maximumNumberResults) * 4;
    bool globalCapReached = false;
    while (coresWound < expectedWoundCores && whileIteration < maxWhileIterations && evaluatedCores.size() < maxEvaluatedCores && !globalCapReached) {
        whileIteration++;
        requestedCores += 20;  // Linear growth instead of exponential
        auto masMagneticsWithCore = coreAdviser.get_advised_core(inputs, coreWeights, requestedCores);

        if (previouslyObtainedCores == masMagneticsWithCore.size()) {
            break;
        }
        previouslyObtainedCores = masMagneticsWithCore.size();

        for (auto& [mas, coreScoring] : masMagneticsWithCore) {
            auto coreNameOpt = mas.get_magnetic().get_core().get_name();
            if (!coreNameOpt) {
                continue;
            }
            std::string coreName = coreNameOpt.value();
            
            if (std::find(evaluatedCores.begin(), evaluatedCores.end(), coreName) != evaluatedCores.end()) {
                continue;
            }
            else {
                evaluatedCores.push_back(coreName);
            }

            // Check performance limit
            if (evaluatedCores.size() >= maxEvaluatedCores) {
                logEntry("Reached maxEvaluatedCores limit (" + std::to_string(maxEvaluatedCores) + ")", "MagneticAdviser", 2);
                break;
            }

            logEntry("core: " + coreName, "MagneticAdviser", 2);
            logEntry("Getting coil", "MagneticAdviser", 2);
            std::vector<std::pair<size_t, double>> usedNumberSectionsAndMargin;
            auto masMagneticsWithCoreAndCoil = coilAdviser.get_advised_coil(mas, std::max(2.0, ceil(double(maximumNumberResults) / masMagneticsWithCore.size())));

            if (masMagneticsWithCoreAndCoil.size() > 0) {
                logEntry("Core wound!", "MagneticAdviser", 2);
                coresWound++;
            }
            size_t processedCoils = 0;
            for (auto mas : masMagneticsWithCoreAndCoil) {

                auto sectionsOpt = mas.get_magnetic().get_coil().get_sections_description();
                if (!sectionsOpt || sectionsOpt->empty()) {
                    continue;
                }
                // A candidate whose turns never wound (sections laid out but
                // turn placement bailed — e.g. a toroid whose turns can't
                // physically fit) must not enter the result pool: downstream
                // simulate()/painter throw COIL_NOT_PROCESSED on the missing
                // turns. The retry-without-toroids loop below already applies
                // this same guard; the main loop was missing it, so an
                // unwindable core leaked through as a turns-less magnetic.
                if (!mas.get_magnetic().get_coil().get_turns_description()) {
                    continue;
                }
                size_t numberSections = sectionsOpt->size();

                double margin = 0;
                if ((*sectionsOpt)[0].get_margin()) {
                    margin = Coil::resolve_margin((*sectionsOpt)[0])[0];
                }
                std::pair<size_t, double> numberSectionsAndMarginCombination = {numberSections, margin};
                if (std::find(usedNumberSectionsAndMargin.begin(), usedNumberSectionsAndMargin.end(), numberSectionsAndMarginCombination) != usedNumberSectionsAndMargin.end()) {
                    continue;
                }

                if (previousCoilIncludeAdditionalCoordinates) {
                    settings.set_coil_include_additional_coordinates(previousCoilIncludeAdditionalCoordinates);
                    mas.get_mutable_magnetic().get_mutable_coil().delimit_and_compact();
                    settings.set_coil_include_additional_coordinates(false);
                }
                try {
                    mas = magneticSimulator.simulate(mas);
                } catch (const std::exception& e) {
                    logEntry(std::string("MagneticAdviser: skipping candidate, simulate failed: ") + e.what(), "MagneticAdviser", 2);
                    continue;
                }

                processedCoils++;

                // Final saturation gate on the ASSEMBLED magnetic. score_magnetics
                // only scores (it discards the validity flag), and the CoreAdviser
                // saturation gate ran on the SEED turns — but the coil adviser /
                // loss optimisation can finalise a higher turn count, lowering the
                // gap-aware saturation current below the margin. Drop inductor
                // designs whose FINAL isat no longer clears margin * ipeak (the same
                // identity downstream realism checks use). Transformers are excluded:
                // their flux is voltage-driven and more turns LOWERS B.
                {
                    bool isInductor = is_inductor(mas.get_mutable_inputs());
                    if (isInductor) {
                        double saturationMargin = settings.get_core_adviser_saturation_margin();
                        bool saturates = false;
                        for (auto& op : mas.get_mutable_inputs().get_operating_points()) {
                            auto excitation = op.get_excitations_per_winding()[0];
                            if (!excitation.get_current() || !excitation.get_current()->get_processed()
                                || !excitation.get_current()->get_processed()->get_peak()) {
                                continue;
                            }
                            double currentPeak = excitation.get_current()->get_processed()->get_peak().value();
                            double temperature = op.get_conditions().get_ambient_temperature();
                            double saturationCurrent = mas.get_mutable_magnetic().calculate_saturation_current(temperature);
                            if (saturationCurrent < saturationMargin * currentPeak) {
                                saturates = true;
                                break;
                            }
                        }
                        if (saturates) {
                            logEntry("MagneticAdviser: dropping '" + mas.get_mutable_magnetic().get_reference()
                                     + "' — final saturation current below margin", "MagneticAdviser", 2);
                            continue;
                        }
                    }
                }

                masData.push_back(mas);
                if (masData.size() >= globalCandidateCap) {
                    logEntry("Reached globalCandidateCap (" + std::to_string(globalCandidateCap) + ")", "MagneticAdviser", 2);
                    globalCapReached = true;
                    break;
                }
                if (processedCoils >= perCoreCoilCap) {
                    usedNumberSectionsAndMargin.push_back(numberSectionsAndMarginCombination);
                    break;
                }
            }
            if (globalCapReached) {
                break;
            }
            if (coresWound >= expectedWoundCores) {
                break;
            }
        }
    }

    logEntry("Found " + std::to_string(masData.size()) + " magnetics", "MagneticAdviser", 2);
    
    auto masMagneticsWithScoring = score_magnetics(masData, filterFlow);

        sort(masMagneticsWithScoring.begin(), masMagneticsWithScoring.end(), [](std::pair<Mas, double>& b1, std::pair<Mas, double>& b2) {
            if (b1.second != b2.second) {
                return b1.second > b2.second;
            }
            // Deterministic tiebreaker: lexicographic reference. Without this,
            // parts with equal totalScoring (common when one filter dominates
            // and other filter contributions saturate) get reordered between
            // runs by std::sort's unstable behaviour, breaking ranking
            // determinism in the UI and tests.
            return b1.first.get_mutable_magnetic().get_reference() < b2.first.get_mutable_magnetic().get_reference();
        });

    if (masMagneticsWithScoring.size() > maximumNumberResults) {
        masMagneticsWithScoring = std::vector<std::pair<Mas, double>>(masMagneticsWithScoring.begin(), masMagneticsWithScoring.end() - (masMagneticsWithScoring.size() - maximumNumberResults));
    }

    // Retry without toroids if toroids were enabled but no results found
    if (masMagneticsWithScoring.empty() && toroidsOriginallyEnabled) {
        logEntry("No magnetics found with toroids enabled. Retrying without toroids...", "MagneticAdviser", 0);
        settings.set_use_toroidal_cores(false);
        clear_loaded_cores();
        // Reset evaluated cores to allow re-evaluation
        evaluatedCores.clear();
        coresWound = 0;
        masData.clear();
        
        // Retry the core search
        whileIteration = 0;
        requestedCores = expectedWoundCores;
        previouslyObtainedCores = SIZE_MAX;
        globalCapReached = false;
        while (coresWound < expectedWoundCores && whileIteration < maxWhileIterations && evaluatedCores.size() < maxEvaluatedCores && !globalCapReached) {
            whileIteration++;
            requestedCores += 20;
            auto masMagneticsWithCore = coreAdviser.get_advised_core(inputs, coreWeights, requestedCores);

            if (previouslyObtainedCores == masMagneticsWithCore.size()) {
                break;
            }
            previouslyObtainedCores = masMagneticsWithCore.size();
            
            for (auto& [mas, coreScoring] : masMagneticsWithCore) {
                auto coreNameOpt = mas.get_magnetic().get_core().get_name();
                if (!coreNameOpt) {
                    continue;
                }
                std::string coreName = coreNameOpt.value();
                
                if (std::find(evaluatedCores.begin(), evaluatedCores.end(), coreName) != evaluatedCores.end()) {
                    continue;
                }
                else {
                    evaluatedCores.push_back(coreName);
                }

                if (evaluatedCores.size() >= maxEvaluatedCores) {
                    logEntry("Reached maxEvaluatedCores limit (" + std::to_string(maxEvaluatedCores) + ")", "MagneticAdviser", 2);
                    break;
                }

                logEntry("core: " + coreName, "MagneticAdviser", 2);
                logEntry("Getting coil", "MagneticAdviser", 2);
                std::vector<std::pair<size_t, double>> usedNumberSectionsAndMargin;
                auto masMagneticsWithCoreAndCoil = coilAdviser.get_advised_coil(mas, std::max(2.0, ceil(double(maximumNumberResults) / masMagneticsWithCore.size())));
                if (masMagneticsWithCoreAndCoil.size() > 0) {
                    logEntry("Core wound!", "MagneticAdviser", 2);
                    coresWound++;
                }
                size_t processedCoils = 0;

                for (auto& masWithCoil : masMagneticsWithCoreAndCoil) {
                    if (masWithCoil.get_magnetic().get_coil().get_turns_description()) {
                        masData.push_back(masWithCoil);
                        processedCoils++;
                        if (masData.size() >= globalCandidateCap) {
                            logEntry("Reached globalCandidateCap (" + std::to_string(globalCandidateCap) + ") in retry", "MagneticAdviser", 2);
                            globalCapReached = true;
                            break;
                        }
                        if (processedCoils >= perCoreCoilCap) {
                            break;
                        }
                    }
                }
            }
            if (globalCapReached) {
                break;
            }
            if (coresWound >= expectedWoundCores) {
                break;
            }
        }
        
        logEntry("Found " + std::to_string(masData.size()) + " magnetics without toroids", "MagneticAdviser", 2);
        masMagneticsWithScoring = score_magnetics(masData, filterFlow);
        
        sort(masMagneticsWithScoring.begin(), masMagneticsWithScoring.end(), [](std::pair<Mas, double>& b1, std::pair<Mas, double>& b2) {
            return b1.second > b2.second;
        });

        if (masMagneticsWithScoring.size() > maximumNumberResults) {
            masMagneticsWithScoring = std::vector<std::pair<Mas, double>>(masMagneticsWithScoring.begin(), masMagneticsWithScoring.end() - (masMagneticsWithScoring.size() - maximumNumberResults));
        }
    }

    // coil_include_additional_coordinates and the three core-filter settings
    // (plus the cached core DB) are restored on scope exit by
    // coilIncludeCoordinatesGuard and coreFilterSettingsRestorer (declared
    // above), so they survive exceptions thrown anywhere in the body.
    return masMagneticsWithScoring;
}

std::vector<std::pair<Mas, double>> MagneticAdviser::get_advised_magnetic(Inputs inputs, std::vector<Magnetic> catalogueMagnetics, size_t maximumNumberResults, bool strict) {
    return get_advised_magnetic(inputs, catalogueMagnetics, _defaultCatalogueMagneticFilterFlow, maximumNumberResults, strict);
}

std::vector<std::pair<Mas, double>> MagneticAdviser::get_advised_magnetic(Inputs inputs, std::vector<Magnetic> catalogueMagnetics, std::vector<MagneticFilterOperation> filterFlow, size_t maximumNumberResults, bool strict) {
    std::vector<Mas> catalogueMagneticsWithInputs;
    for (auto magnetic : catalogueMagnetics) {
        if (inputs.get_operating_points().size() > 0 && magnetic.get_mutable_coil().get_functional_description().size() != inputs.get_operating_points()[0].get_excitations_per_winding().size()) {
            continue;
        }
        Mas mas;
        mas.set_inputs(inputs);
        mas.set_magnetic(magnetic);
        catalogueMagneticsWithInputs.push_back(mas);
    }

    return get_advised_magnetic(catalogueMagneticsWithInputs, filterFlow, maximumNumberResults, strict);
}

std::vector<std::pair<Mas, double>> MagneticAdviser::get_advised_magnetic(Inputs inputs, const std::map<std::string, Magnetic>& catalogueMagnetics, std::vector<MagneticFilterOperation> filterFlow, size_t maximumNumberResults, bool strict) {
    std::vector<Mas> catalogueMagneticsWithInputs;
    catalogueMagneticsWithInputs.reserve(catalogueMagnetics.size());
    for (const auto& [reference, magnetic] : catalogueMagnetics) {
        if (inputs.get_operating_points().size() > 0 && magnetic.get_coil().get_functional_description().size() != inputs.get_operating_points()[0].get_excitations_per_winding().size()) {
            continue;
        }
        Mas mas;
        mas.set_inputs(inputs);
        mas.set_magnetic(magnetic);
        catalogueMagneticsWithInputs.push_back(mas);
    }

    return get_advised_magnetic(catalogueMagneticsWithInputs, filterFlow, maximumNumberResults, strict);
}

std::vector<std::pair<Mas, double>> MagneticAdviser::get_advised_magnetic(std::vector<Mas> catalogueMagneticsWithInputs, std::vector<MagneticFilterOperation> filterFlow, size_t maximumNumberResults, bool strict) {

    // No candidate magnetics survived upstream filtering — e.g. the requested
    // winding count has no matching parts in the catalogue (the callers above
    // skip any magnetic whose functionalDescription size differs from the
    // operating point's excitations-per-winding count). With an empty vector the
    // `catalogueMagneticsWithInputs[0]` access below reads past the end and
    // segfaults in the Inputs copy constructor. "No matching candidate" is a
    // legitimate outcome, so return an empty result set instead of crashing.
    if (catalogueMagneticsWithInputs.empty()) {
        return {};
    }

    // Clear the process-global _scorings map so stale entries from previous adviser
    // calls (same process, same catalog) do not pollute this run's min/max
    // normalization. Without this, top scores never reach 1.0 and identical inputs
    // can produce different rankings between runs.
    clear_scoring();

    load_filter_flow(filterFlow, catalogueMagneticsWithInputs[0].get_inputs());
    std::vector<MagneticFilterOperation> strictlyRequiredFilterFlow;
    std::vector<MagneticFilterOperation> nonStrictlyRequiredFilterFlow;
    std::vector<Mas> validMas;
    MagneticSimulator magneticSimulator;

    std::map<std::string, double> scoringPerReference;
    std::vector<Mas> catalogueMasWithStriclyRequirementsPassed;
    for (auto filterConfiguration : filterFlow) {
        if (filterConfiguration.get_strictly_required()) {
            strictlyRequiredFilterFlow.push_back(filterConfiguration);
            // break;
        }
        else {
            nonStrictlyRequiredFilterFlow.push_back(filterConfiguration);
        }
    }
    // Detect data-corruption-class failures vs per-core failures.
    // A throw is "per-core" when it depends on the specific core
    // (geometry mismatch, bobbin lookup miss, etc.) — different
    // candidates raise different exceptions. A throw is
    // "data-class" when every candidate hits the SAME message
    // (typically an operating-point malformation: zero harmonics,
    // missing excitation, etc.); retrying every candidate is
    // useless and just burns wall-clock to no end. Bail out after
    // a small streak of identical messages so callers see a
    // diagnosable error instead of timing out.
    static constexpr size_t kIdenticalThrowBailThreshold = 8;
    std::string previousThrowMessage;
    size_t identicalThrowStreak = 0;
    for (size_t index = 0; index < catalogueMagneticsWithInputs.size(); ++index) {
        auto mas = catalogueMagneticsWithInputs[index];
        std::vector<Outputs> outputs;
        auto inputs = mas.get_inputs();
        auto magnetic = mas.get_magnetic();
        bool validMagnetic = true;
        for (auto filterConfiguration : strictlyRequiredFilterFlow) {
            MagneticFilters filterEnum = filterConfiguration.get_filter();

            try {
                auto [valid, scoring] = _filters[filterEnum]->evaluate_magnetic(&magnetic, &inputs, &outputs);
                add_scoring(magnetic.get_reference(), filterEnum, scoring);
                if (strict) {
                    validMagnetic &= valid;
                    if (!valid) {
                        break;
                    }
                }
                // Successful evaluation resets the streak.
                identicalThrowStreak = 0;
                previousThrowMessage.clear();
            }
            catch (const std::exception& e) {
                logEntry(std::string("MagneticAdviser: strict filter ") + std::string(magic_enum::enum_name(filterEnum)) + " threw, rejecting magnetic: " + e.what(), "MagneticAdviser", 2);
                std::string thisMsg = e.what();
                if (thisMsg == previousThrowMessage) {
                    ++identicalThrowStreak;
                } else {
                    previousThrowMessage = thisMsg;
                    identicalThrowStreak = 1;
                }
                validMagnetic = false;
                break;
            }
        }
        if (identicalThrowStreak >= kIdenticalThrowBailThreshold) {
            throw InvalidInputException(
                ErrorCode::INVALID_INPUT,
                std::string("MagneticAdviser: aborted after ") +
                std::to_string(identicalThrowStreak) +
                " consecutive candidates failed with identical error '" +
                previousThrowMessage +
                "'. The failure is not core-specific — the operating point or " +
                "inputs are malformed (e.g. zero-harmonic excitation). Fix " +
                "process_design_requirements / process_operating_points for " +
                "this topology rather than continuing the candidate search."
            );
        }

        if (validMagnetic) {
            Mas resultMas;
            resultMas.set_magnetic(magnetic);
            resultMas.set_inputs(inputs);
            resultMas.set_outputs(outputs);
            catalogueMasWithStriclyRequirementsPassed.push_back(resultMas);
        }
    }

    for (size_t index = 0; index < catalogueMasWithStriclyRequirementsPassed.size(); ++index) {
        auto mas = catalogueMasWithStriclyRequirementsPassed[index];
        // Seed `outputs` from the strict-loop's accumulated outputs so the
        // final resultMas carries them. Previously this re-initialised to an
        // empty vector and the strict-loop outputs (e.g. core/winding losses
        // from LOSSES_NO_PROXIMITY) were discarded.
        auto outputs = mas.get_outputs();
        auto inputs = mas.get_inputs();
        auto magnetic = mas.get_magnetic();
        bool valid = true;
        // Iterate only the non-strict filters here. The strict ones were already
        // evaluated and scored in the previous loop; re-running them is wasted work
        // (and ~2x slowdown for filter flows dominated by strict filters such as CMC).
        // Re-running is otherwise harmless because add_scoring overwrites by key.
        for (auto filterConfiguration : nonStrictlyRequiredFilterFlow) {
            MagneticFilters filterEnum = filterConfiguration.get_filter();

            try {
                // Loop B is intentionally score-only: the per-filter `valid` flag is
                // discarded for non-strict filters. We rely on the outer `valid`
                // staying true so this magnetic survives into `validMas`.
                //
                // Do NOT change this to `valid &= filterValid` without also handling
                // the case where every part is rejected: the empty-validMas branch
                // below recurses with `strict=false`, and a non-strict filter that
                // still rejects everything would recurse forever and stack-overflow.
                auto [filterValid, scoring] = _filters[filterEnum]->evaluate_magnetic(&magnetic, &inputs, &outputs);
                (void)filterValid;
                add_scoring(magnetic.get_reference(), filterEnum, scoring);
            }
            catch (const std::exception& e) {
                logEntry(std::string("MagneticAdviser: non-strict filter ") + std::string(magic_enum::enum_name(filterEnum)) + " threw, rejecting magnetic: " + e.what(), "MagneticAdviser", 2);
                valid = false;
                break;
            }
        }

        if (valid) {
            Mas resultMas;
            resultMas.set_magnetic(magnetic);
            resultMas.set_inputs(inputs);
            resultMas.set_outputs(outputs);
            validMas.push_back(resultMas);
        }
    }

    auto scoringsPerReferencePerFilter = get_scorings();

    // Build weight lookup from filter flow for weighted scoring
    std::map<MagneticFilters, double> filterWeights;
    double totalWeight = 0;
    for (auto& filterConfiguration : filterFlow) {
        filterWeights[filterConfiguration.get_filter()] = filterConfiguration.get_weight();
        totalWeight += filterConfiguration.get_weight();
    }

    std::vector<std::pair<Mas, double>> masMagneticsWithScoring;

    if (validMas.size() > 0) {

        for (auto mas : validMas) {
            auto reference = mas.get_mutable_magnetic().get_reference();
            double totalScoring = 0;
            double usedWeight = 0;
            for (auto [filter, scoring] : scoringsPerReferencePerFilter[reference]) {
                double weight = filterWeights.count(filter) ? filterWeights[filter] : 1.0;
                totalScoring += scoring * weight;
                usedWeight += weight;
            }
            // Normalize by total weight used (weighted average)
            if (usedWeight > 0) {
                totalScoring /= usedWeight;
            }
            masMagneticsWithScoring.push_back({mas, totalScoring});
        }

        sort(masMagneticsWithScoring.begin(), masMagneticsWithScoring.end(), [](std::pair<Mas, double>& b1, std::pair<Mas, double>& b2) {
            return b1.second > b2.second;
        });
        if (masMagneticsWithScoring.size() > maximumNumberResults) {
            masMagneticsWithScoring = std::vector<std::pair<Mas, double>>(masMagneticsWithScoring.begin(), masMagneticsWithScoring.end() - (masMagneticsWithScoring.size() - maximumNumberResults));
        }

        if (_simulateResults) {
            std::vector<std::pair<Mas, double>> masMagneticsWithScoringSimulated;
            for (auto [mas, scoring] : masMagneticsWithScoring) {
                try {
                    mas = magneticSimulator.simulate(mas, true);
                } catch (const std::exception& e) {
                    logEntry(std::string("MagneticAdviser: skipping final-simulate candidate: ") + e.what(), "MagneticAdviser", 2);
                    continue;
                }
                masMagneticsWithScoringSimulated.push_back({mas, scoring});
            }

            return masMagneticsWithScoringSimulated;
        }
        else {
            return masMagneticsWithScoring;
        }

    }
    else {
        if (catalogueMasWithStriclyRequirementsPassed.size() > 0) {
            return get_advised_magnetic(catalogueMasWithStriclyRequirementsPassed, filterFlow, maximumNumberResults, false);
        }
        return {};
    }
}

std::vector<std::pair<Mas, double>> MagneticAdviser::score_magnetics(std::vector<Mas> masMagnetics, std::vector<MagneticFilterOperation> filterFlow) {
    std::vector<std::pair<Mas, double>> masMagneticsWithScoring;

    // Return early if no magnetics to score
    if (masMagnetics.empty()) {
        return masMagneticsWithScoring;
    }

    for (auto mas : masMagnetics) {
        masMagneticsWithScoring.push_back({mas, 0.0});
    }
    for (auto filterConfiguration : filterFlow) {
        MagneticFilters filterEnum = filterConfiguration.get_filter();

        // Check if filter exists in the map before accessing
        auto filterIt = _filters.find(filterEnum);
        if (filterIt == _filters.end() || !filterIt->second) {
            continue;  // Skip filters that aren't loaded
        }

        std::vector<double> scorings;
        for (auto mas : masMagnetics) {
            auto [valid, scoring] = filterIt->second->evaluate_magnetic(&mas.get_mutable_magnetic(), &mas.get_mutable_inputs());
            scorings.push_back(scoring);
            add_scoring(mas.get_mutable_magnetic().get_reference(), filterEnum, scoring);
        }
        if (masMagneticsWithScoring.size() > 0) {
            normalize_scoring(&masMagneticsWithScoring, scorings, filterConfiguration);
        }
    }
    
    // OPTIMIZATION: Apply stack penalty - prioritize single stack designs
    // Penalty increases with number of stacks: 1 stack = 0%, 2 stacks = 15%, 3 stacks = 30%, etc.
    for (auto& [mas, scoring] : masMagneticsWithScoring) {
        auto numStacksOpt = mas.get_magnetic().get_core().get_functional_description().get_number_stacks();
        if (numStacksOpt && numStacksOpt.value() > 1) {
            double stackPenalty = (numStacksOpt.value() - 1) * 0.15;  // 15% penalty per additional stack
            scoring *= (1.0 - stackPenalty);
        }
    }
    
    return masMagneticsWithScoring;
}


void MagneticAdviser::preview_magnetic(Mas mas) {
    std::string text = "";

    text += "Core shape: " + mas.get_mutable_magnetic().get_mutable_core().get_shape_name() + "\n";
    text += "Core material: " + mas.get_mutable_magnetic().get_mutable_core().get_material_name() + "\n";
    if (mas.get_mutable_magnetic().get_mutable_core().get_functional_description().get_gapping().size() > 0) {
        text += "Core gap: " + std::to_string(mas.get_mutable_magnetic().get_mutable_core().get_functional_description().get_gapping()[0].get_length()) + "\n";
    }
    auto numberStacksOpt = mas.get_mutable_magnetic().get_mutable_core().get_functional_description().get_number_stacks();
    if (numberStacksOpt) {
        text += "Core stacks: " + std::to_string(numberStacksOpt.value()) + "\n";
    }
    for (size_t windingIndex = 0; windingIndex < mas.get_mutable_magnetic().get_mutable_coil().get_functional_description().size(); ++windingIndex) {
        auto winding = mas.get_mutable_magnetic().get_mutable_coil().get_functional_description()[windingIndex];
        auto wire = OpenMagnetics::Coil::resolve_wire(winding);
        text += "Winding: " + winding.get_name() + "\n";
        text += "\tNumber Turns: " + std::to_string(winding.get_number_turns()) + "\n";
        text += "\tNumber Parallels: " + std::to_string(winding.get_number_parallels()) + "\n";
        text += "\tWire: " + to_string(wire.get_type());
        if (wire.get_standard()) {
            text += " " + to_string(wire.get_standard().value());
        }
        if (wire.get_name()) {
            text += " " + wire.get_name().value();
        }
        text += "\n";
    }

    for (size_t operatingPointIndex = 0; operatingPointIndex < mas.get_outputs().size(); ++operatingPointIndex) {
        auto output = mas.get_outputs()[operatingPointIndex];
        text += "Operating Point: " + std::to_string(operatingPointIndex + 1) + "\n";
        text += "\tMagnetizing Inductance: " + std::to_string(resolve_dimensional_values(output.get_inductance()->get_magnetizing_inductance().get_magnetizing_inductance())) + "\n";
        
        if (output.get_core_losses()) {
            auto coreLosses = output.get_core_losses().value();
            text += "\tCore losses: " + std::to_string(coreLosses.get_core_losses()) + "\n";
            
            auto fluxDensity = coreLosses.get_magnetic_flux_density();
            if (fluxDensity && fluxDensity->get_processed() && fluxDensity->get_processed()->get_peak()) {
                text += "\tMagnetic flux density: " + std::to_string(fluxDensity->get_processed()->get_peak().value()) + "\n";
            }
            
            if (coreLosses.get_temperature() > 0) {
                text += "\tCore temperature: " + std::to_string(coreLosses.get_temperature()) + "\n";
            }
        }
        
        if (output.get_winding_losses()) {
            auto windingLosses = output.get_winding_losses().value();
            text += "\tWinding losses: " + std::to_string(windingLosses.get_winding_losses()) + "\n";
            
            if (windingLosses.get_winding_losses_per_winding()) {
                auto lossesPerWinding = windingLosses.get_winding_losses_per_winding().value();
                for (size_t windingIndex = 0; windingIndex < lossesPerWinding.size(); ++windingIndex) {
                    auto winding = mas.get_mutable_magnetic().get_mutable_coil().get_functional_description()[windingIndex];
                    auto windingLoss = lossesPerWinding[windingIndex];
                    text += "\t\tLosses for winding: " + winding.get_name() + "\n";
                    double skinEffectLosses = 0;
                    double proximityEffectLosses = 0;
                    
                    if (windingLoss.get_skin_effect_losses()) {
                        auto skinLosses = windingLoss.get_skin_effect_losses().value();
                        for (size_t i = 0; i < skinLosses.get_losses_per_harmonic().size(); ++i) {
                            skinEffectLosses += skinLosses.get_losses_per_harmonic()[i];
                        }
                    }
                    
                    if (windingLoss.get_proximity_effect_losses()) {
                        auto proximityLosses = windingLoss.get_proximity_effect_losses().value();
                        for (size_t i = 0; i < proximityLosses.get_losses_per_harmonic().size(); ++i) {
                            proximityEffectLosses += proximityLosses.get_losses_per_harmonic()[i];
                        }
                    }

                    if (windingLosses.get_dc_resistance_per_winding() && windingIndex < windingLosses.get_dc_resistance_per_winding().value().size()) {
                        text += "\t\t\tDC resistance: " + std::to_string(windingLosses.get_dc_resistance_per_winding().value()[windingIndex]) + "\n";
                    }
                    
                    if (windingLoss.get_ohmic_losses()) {
                        text += "\t\t\tOhmic losses: " + std::to_string(windingLoss.get_ohmic_losses().value().get_losses()) + "\n";
                    }
                    
                    text += "\t\t\tSkin effect losses: " + std::to_string(skinEffectLosses) + "\n";
                    text += "\t\t\tProximity effect losses: " + std::to_string(proximityEffectLosses) + "\n";

                    if (windingIndex > 0) {
                        auto leakageInductance = output.get_inductance()->get_leakage_inductance();
                        if (leakageInductance && windingIndex - 1 < leakageInductance->get_leakage_inductance_per_winding().size()) {
                            auto nominalOpt = leakageInductance->get_leakage_inductance_per_winding()[windingIndex - 1].get_nominal();
                            if (nominalOpt) {
                                double value = nominalOpt.value();
                                text += "\t\t\tLeakage inductance referred to primary: " + std::to_string(value) + "\n";
                            }
                        }
                    }
                }
            }
        }
    }
}


std::map<std::string, std::map<MagneticFilters, double>> MagneticAdviser::get_scorings(){
    std::map<std::string, std::map<MagneticFilters, double>> swappedScorings;
    for (auto& [filter, aux] : _scorings) {

        MagneticFilterOperation magneticFilterOperation;
        for (auto loadedFilter : _loadedFilterFlow) {
            if (loadedFilter.get_filter() == filter) {
                magneticFilterOperation = loadedFilter;
                break;
            }
        }

        auto normalizedScorings = OpenMagnetics::normalize_scoring(aux, magneticFilterOperation);
        for (auto& [name, scoring] : aux) {
            swappedScorings[name][filter] = normalizedScorings[name];
        }
    }
    return swappedScorings;
}

} // namespace OpenMagnetics

bool is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

// int main(int argc, char* argv[]) {
//     if (argc == 1) {
//         throw std::invalid_argument("Missing inputs file");
//     }
//     else {
//         int numberMagnetics = 1;
//         std::filesystem::path inputFilepath = argv[1];
//         auto outputFilePath = std::filesystem::path {std::source_location::current().file_name()}.parent_path().append("..").append("output");
//         if (argc >= 3) {
//             if (is_number(argv[2])) {
//                 numberMagnetics = std::stoi(argv[2]);
//             }
//             else {
//                 outputFilePath = argv[2];
//             }
//         }
//         if (argc >= 4) {
//             if (!is_number(argv[2]) && is_number(argv[3])) {
//                 numberMagnetics = std::stoi(argv[3]);
//             }
//         }

//         std::ifstream f(inputFilepath);
//         json masJson;
//         std::string str;
//         if(static_cast<bool>(f)) {
//             std::ostringstream ss;
//             ss << f.rdbuf(); // reading data
//             masJson = json::parse(ss.str());
//         }

//         if (masJson["inputs"]["designRequirements"]["insulation"]["cti"] == "Group IIIa") {
//             masJson["inputs"]["designRequirements"]["insulation"]["cti"] = "Group IIIA";
//         }
//         if (masJson["inputs"]["designRequirements"]["insulation"]["cti"] == "Group IIIb") {
//             masJson["inputs"]["designRequirements"]["insulation"]["cti"] = "Group IIIB";
//         }

//         OpenMagnetics::Inputs inputs(masJson["inputs"]);

//         auto outputFolder = inputFilepath.parent_path();


//         OpenMagnetics::MagneticAdviser MagneticAdviser;
//         auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, numberMagnetics);
//         for (size_t i = 0; i < masMagnetics.size(); ++i){
//             OpenMagnetics::MagneticAdviser::preview_magnetic(masMagnetics[i].first);

//             std::filesystem::path outputFilename = outputFilePath;
//             outputFilename += inputFilepath.filename();
//             outputFilename += "_design_" + std::to_string(i) + ".json";
//             std::ofstream outputFile(outputFilename);
            
//             json output;
//             to_json(output, masMagnetics[i].first);
//             outputFile << output;
            
//             outputFile.close();

//             outputFilename.replace_extension("svg");
//             OpenMagnetics::Painter painter(outputFilename);
//             settings.set_painter_mode(OpenMagnetics::PainterModes::CONTOUR);

//             settings.set_painter_number_points_x(20);
//             settings.set_painter_number_points_y(20);
//             settings.set_painter_include_fringing(false);
//             settings.set_painter_mirroring_dimension(0);

//             painter.paint_magnetic_field(masMagnetics[i].first.get_mutable_inputs().get_operating_point(0), masMagnetics[i].first.get_mutable_magnetic());
//             painter.paint_core(masMagnetics[i].first.get_mutable_magnetic());
//             painter.paint_bobbin(masMagnetics[i].first.get_mutable_magnetic());
//             painter.paint_coil_turns(masMagnetics[i].first.get_mutable_magnetic());
//             painter.export_svg();
//         }


//     }
// }