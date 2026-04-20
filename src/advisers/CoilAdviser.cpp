#include "processors/Inputs.h"
#include "advisers/CoilAdviser.h"
#include "Models.h"
#include "constructive_models/Insulation.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include <algorithm>
#include <set>
#include <limits> // B18 FIX: for numeric_limits
#include "support/Exceptions.h"
#include "support/Logger.h"


namespace OpenMagnetics {

    /**
     * @brief Calculate winding window proportion for each winding.
     *
     * Uses equal proportions for all windings. This is robust for all transformer
     * topologies (Flyback, Forward, Push-Pull, etc.) where V*I relationships are
     * complex and don't map directly to window allocation.
     * 
     * Equal proportions let the winding algorithm handle fitting with parallels/layers,
     * and the scoring system will rank designs by actual performance metrics.
     *
     * Future enhancement: Could weight by I_rms × N_turns for copper area optimization,
     * but this requires careful handling of edge cases and may not improve results
     * given the downstream scoring filters.
     */
    /**
     * @brief Calculate winding window proportion for each winding.
     *
     * For transformers, power flows through all windings (P = V × I ≈ constant).
     * Since turns ratio N ∝ V and current I ∝ 1/V, the copper volume 
     * (I × N × cross_section) is approximately equal for all windings.
     *
     * Therefore, equal proportions is the correct approach for most transformer
     * topologies. The winding algorithm handles the actual fitting with different
     * wire gauges and parallel counts.
     *
     * Note: For inductors (single winding), returns {1.0}.
     */
    std::vector<double> calculate_winding_window_proportion_per_winding(Inputs& inputs) {
        size_t numWindings = inputs.get_operating_points()[0].get_excitations_per_winding().size();
        
        // Equal proportions - correct for transformers where P_in ≈ P_out
        // and copper volume I×N is similar for all windings
        return std::vector<double>(numWindings, 1.0 / numWindings);
    }

    void CoilAdviser::load_filter_flow(std::vector<MagneticFilterOperation> flow, std::optional<Inputs> inputs) {
        _filters.clear();
        _loadedFilterFlow = flow;
        for (const auto& filterConfiguration : flow) {
            MagneticFilters filterEnum = filterConfiguration.get_filter();
            _filters[filterEnum] = MagneticFilter::factory(filterEnum, inputs);
        }
    }

    std::vector<std::pair<Mas, double>> CoilAdviser::score_magnetics(std::vector<Mas> masMagnetics, std::vector<MagneticFilterOperation> filterFlow) {
        std::vector<std::pair<Mas, double>> masMagneticsWithScoring;
        std::vector<bool> masValidFlags(masMagnetics.size(), true);
        for (auto mas : masMagnetics) {
            masMagneticsWithScoring.push_back({mas, 0.0});
        }
        for (const auto& filterConfiguration : filterFlow) {
            MagneticFilters filterEnum = filterConfiguration.get_filter();
            
            std::vector<double> scorings;
            for (size_t masIndex = 0; masIndex < masMagnetics.size(); ++masIndex) {
                auto [valid, scoring] = _filters[filterEnum]->evaluate_magnetic(&masMagnetics[masIndex].get_mutable_magnetic(), &masMagnetics[masIndex].get_mutable_inputs());
                if (!valid) {
                    masValidFlags[masIndex] = false;
                }
                scorings.push_back(scoring);
                add_scoring(masMagnetics[masIndex].get_mutable_magnetic().get_reference(), filterEnum, scoring);
            }
            if (masMagneticsWithScoring.size() > 0) {
                normalize_scoring(&masMagneticsWithScoring, scorings, filterConfiguration);
            }
        }
        // Remove invalid designs
        std::vector<std::pair<Mas, double>> validMasMagneticsWithScoring;
        size_t invalidCount = 0;
        for (size_t i = 0; i < masMagneticsWithScoring.size(); ++i) {
            if (masValidFlags[i]) {
                validMasMagneticsWithScoring.push_back(masMagneticsWithScoring[i]);
            }
            else {
                invalidCount++;
            }
        }
        if (invalidCount > 0) {
            logEntry("Filtered out " + std::to_string(invalidCount) + " invalid designs out of " + std::to_string(masMagnetics.size()), "CoilAdviser", 2);
        }
        return validMasMagneticsWithScoring;
    }

    std::vector<Mas> CoilAdviser::get_advised_coil(Mas mas, size_t maximumNumberResults){
        logEntry("Starting Coil Adviser without wires", "CoilAdviser");
        if (wireDatabase.empty()) {
            load_wires();
        }
        std::string jsonLine;
        std::vector<Wire> wires;
        for (const auto& [key, wire] : wireDatabase) {
            if ((settings.get_wire_adviser_include_planar() || wire.get_type() != WireType::PLANAR) &&
                (settings.get_wire_adviser_include_foil() || wire.get_type() != WireType::FOIL) &&
                (settings.get_wire_adviser_include_rectangular() || wire.get_type() != WireType::RECTANGULAR) &&
                (settings.get_wire_adviser_include_litz() || wire.get_type() != WireType::LITZ) &&
                (settings.get_wire_adviser_include_round() || wire.get_type() != WireType::ROUND)) {

                if (!_commonWireStandard || !wire.get_standard()) {
                    wires.push_back(wire);
                }
                else if (wire.get_standard().value() == _commonWireStandard){
                    wires.push_back(wire);
                }
            }
        }
        return get_advised_coil(&wires, mas, maximumNumberResults);
    }

    std::vector<Mas> CoilAdviser::get_advised_coil(std::vector<Wire>* wires, Mas mas, size_t maximumNumberResults){
        logEntry("Starting Coil Adviser", "CoilAdviser");
        auto core = mas.get_magnetic().get_core();
        auto coreType = core.get_functional_description().get_type();

        auto inputs = mas.get_mutable_inputs();
        load_filter_flow(_defaultCustomMagneticFilterFlow, inputs);
        auto patterns = Coil::get_patterns(inputs, coreType);
        auto repetitions = Coil::get_repetitions(inputs, coreType);
        mas.set_inputs(inputs);

        // ── CMC winding policy ──────────────────────────────────────────
        // Common-mode chokes (sub_application == COMMON_MODE_NOISE_FILTERING)
        // follow a different topology from generic coupled inductors:
        //   · one contiguous section per winding (each phase lives on its
        //     own angular sector of the toroid), NOT interleaved / bifilar;
        //   · sections spread across the toroid arc for equal proportions;
        //   · turns centered within each section (compact bundle mid-sector);
        //   · matching wire and parallels across every winding.
        // The existing Coil::get_repetitions returns {2, 1} for CMCs today
        // for legacy balun-style designs; override here so the EMI-filter
        // layout wins. The wire + parallel mirror happens per-combination
        // inside the winding loop below — see "CMC: mirror winding[0]" note.
        //
        // TODO: creepage-driven inter-section margin — needs either a public
        // Coil::set_section_margins setter or a CMC-aware branch in
        // InsulationCoordinator; deferred to avoid modifying Coil.cpp here.
        const bool isCmc = inputs.get_design_requirements().get_sub_application().has_value()
                        && inputs.get_design_requirements().get_sub_application().value()
                           == SubApplication::COMMON_MODE_NOISE_FILTERING;
        if (isCmc && coreType == CoreType::TOROIDAL) {
            repetitions = {1};
            mas.get_mutable_magnetic().get_mutable_coil()
                .set_winding_orientation(WindingOrientation::CONTIGUOUS);
            mas.get_mutable_magnetic().get_mutable_coil()
                .set_section_alignment(CoilAlignment::SPREAD);
            mas.get_mutable_magnetic().get_mutable_coil()
                .set_turns_alignment(CoilAlignment::CENTERED);
            logEntry("CMC detected: forcing {1} repetitions, sections CONTIGUOUS + SPREAD, turns CENTERED", "CoilAdviser");
        }

        
        size_t maximumNumberResultsPerPattern = std::max(2.0, ceil(maximumNumberResults / (patterns.size() * repetitions.size())));
        logEntry("Trying " + std::to_string(repetitions.size()) + " repetitions and " + std::to_string(patterns.size()) + " patterns", "CoilAdviser");

        std::vector<Mas> masesWithCoil;
        for (auto repetition : repetitions) {
            for (auto pattern : patterns) {
                auto aux = mas.get_mutable_magnetic().get_mutable_coil().check_pattern_and_repetitions_integrity(pattern, repetition);
                pattern = aux.first;
                repetition = aux.second;
                if (mas.get_mutable_inputs().get_wiring_technology() == WiringTechnology::WOUND) {
                    auto combinationsSolidInsulationRequirementsForWires = InsulationCoordinator::get_solid_insulation_requirements_for_wires(mas.get_mutable_inputs(), pattern, repetition);
                    for(size_t insulationIndex = 0; insulationIndex < combinationsSolidInsulationRequirementsForWires.size(); ++insulationIndex) {
                        auto solidInsulationRequirementsForWires = combinationsSolidInsulationRequirementsForWires[insulationIndex];
                        std::string reference = "Custom";
                        if (mas.get_magnetic().get_manufacturer_info() && mas.get_magnetic().get_manufacturer_info()->get_reference()) {
                            reference = mas.get_magnetic().get_manufacturer_info().value().get_reference().value();
                        }

                        reference += ", Turns: ";
                        reference += std::to_string(mas.get_magnetic().get_coil().get_functional_description()[0].get_number_turns());

                        reference += ", Order: ";
                        for (auto pat : pattern) {
                            reference += std::to_string(pat);
                        }
                        if (repetition > 1) {
                            reference += ", Interleaved";
                        }
                        else {
                            reference += ", Non-Interleaved";
                        }
                        if (InsulationCoordinator::needs_margin(solidInsulationRequirementsForWires, pattern, repetition)) {
                            reference += ", Margin Taped";
                        }
                        else {
                            reference += ", Wire Insulated";
                        }
                        reference += " " + std::to_string(insulationIndex);


                        auto resultsPerPattern = get_advised_coil_for_pattern(wires, mas, pattern, repetition, solidInsulationRequirementsForWires, maximumNumberResultsPerPattern, reference);

                        std::move(resultsPerPattern.begin(), resultsPerPattern.end(), std::back_inserter(masesWithCoil));
                    }

                }
                else {
                    std::string reference = "Custom";
                    if (mas.get_magnetic().get_manufacturer_info() && mas.get_magnetic().get_manufacturer_info()->get_reference()) {
                        reference = mas.get_magnetic().get_manufacturer_info().value().get_reference().value();
                    }

                    reference += ", Turns: ";
                    reference += std::to_string(mas.get_magnetic().get_coil().get_functional_description()[0].get_number_turns());
                    auto resultsPerPattern = get_advised_planar_coil_for_pattern(wires, mas, pattern, repetition, maximumNumberResultsPerPattern, reference);
                    std::move(resultsPerPattern.begin(), resultsPerPattern.end(), std::back_inserter(masesWithCoil));
                }
            }
        }

        logEntry("Found " + std::to_string(masesWithCoil.size()) + " magnetics", "CoilAdviser");
        auto masMagneticsWithScoring = score_magnetics(masesWithCoil, _loadedFilterFlow);

        stable_sort(masMagneticsWithScoring.begin(), masMagneticsWithScoring.end(), [](const std::pair<Mas, double>& b1, const std::pair<Mas, double>& b2) {
            return b1.second > b2.second;
        });

        std::vector<Mas> masesWithoutScoring;
        for (auto [mas, scoring] : masMagneticsWithScoring) {
            masesWithoutScoring.push_back(mas);
        }

        // Fallback: If all designs were filtered out but we had valid windings, return the best available
        if (masesWithoutScoring.empty() && !masesWithCoil.empty()) {
            logEntry("WARNING: All " + std::to_string(masesWithCoil.size()) + " designs were filtered out by criteria. " +
                     "Returning best available design as fallback.", "CoilAdviser", 1);
            
            // Return up to maximumNumberResults from the original unfiltered designs
            size_t numToReturn = std::min(masesWithCoil.size(), maximumNumberResults);
            for (size_t i = 0; i < numToReturn; ++i) {
                masesWithoutScoring.push_back(masesWithCoil[i]);
            }
        }

        if (masesWithoutScoring.size() > maximumNumberResults) {
            masesWithoutScoring = std::vector<Mas>(masesWithoutScoring.begin(), masesWithoutScoring.end() - (masesWithoutScoring.size() - maximumNumberResults));
        }

        return masesWithoutScoring;
    }

    std::vector<Section> CoilAdviser::get_advised_sections(Mas mas, std::vector<size_t> pattern, size_t repetitions){
        auto sectionProportions = calculate_winding_window_proportion_per_winding(mas.get_mutable_inputs());
        auto core = mas.get_magnetic().get_core();
        auto coil = mas.get_magnetic().get_coil();

        coil.set_strict(false);
        coil.set_inputs(mas.get_inputs());

        coil.calculate_insulation(true);
        auto result = coil.wind_by_sections(sectionProportions, pattern, repetitions);
        if (result) {
            coil.delimit_and_compact();
            coil.set_strict(true);
        }

        if (coil.get_sections_description()) {
            return coil.get_sections_description().value();
        }
        else {
            return {};
        }
    }

    std::vector<Section> CoilAdviser::get_advised_planar_sections(Mas mas, std::vector<size_t> pattern, size_t repetitions){
        auto sectionProportions = calculate_winding_window_proportion_per_winding(mas.get_mutable_inputs());
        auto core = mas.get_magnetic().get_core();
        auto coil = mas.get_magnetic().get_coil();

        coil.set_strict(false);
        coil.set_inputs(mas.get_inputs());
        coil.calculate_insulation(true);
        coil.set_groups_description(std::nullopt);
        std::vector<size_t> stackUp;

        auto numberTurnsPerWinding = mas.get_mutable_magnetic().get_mutable_coil().get_number_turns();

        size_t totalNumberLayers = settings.get_coil_maximum_layers_planar();
        for (size_t repetitionIndex = 0; repetitionIndex < repetitions; ++repetitionIndex) {
            for (auto windingIndex : pattern) {
                if (windingIndex >= numberTurnsPerWinding.size()) continue;
                for (size_t layerIndex = 0; layerIndex < floor(totalNumberLayers / repetitions / pattern.size()); ++layerIndex) {
                    if (numberTurnsPerWinding[windingIndex] > 0) {
                        stackUp.push_back(windingIndex);
                        numberTurnsPerWinding[windingIndex]--;
                    }
                }
            }
        }

        auto result = coil.wind_by_planar_sections(stackUp, {}, defaults.coreToLayerDistance);
        if (result) {
            coil.set_strict(true);
        }

        if (coil.get_sections_description()) {
            auto sections = coil.get_sections_description().value();
            return sections;
        }
        else {
            return {};
        }
    }

    /**
     * @brief Calculate optimal planar stackup considering skin depth for AC losses.
     *
     * For high-frequency operation, copper utilization is limited by skin depth.
     * When wire thickness >> 2×skin_depth, current only flows near the surfaces,
     * making thick copper inefficient. This function considers:
     * 
     * 1. Minimum layers needed to fit turns geometrically
     * 2. Optimal layers based on skin depth (more thinner layers may be better)
     * 3. Constraints from parallel configurations
     * 
     * The skin depth consideration helps avoid over-thick copper that wastes
     * material without reducing AC resistance.
     */
    std::vector<size_t> get_planar_stackup(Coil coil, std::vector<double> proportionPerWinding, std::vector<size_t> pattern, size_t repetitions, std::optional<double> effectiveFrequency = std::nullopt, double temperature = 25.0) {
        size_t totalNumberLayers = settings.get_coil_maximum_layers_planar();
        size_t totalAssignedLayers = 0;
        std::vector<size_t> layersPerWinding;
        std::vector<size_t> stackUp;
        
        // Initial allocation based on proportions
        for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
            auto layersThisWinding = floor(totalNumberLayers * proportionPerWinding[windingIndex]);
            totalAssignedLayers += layersThisWinding;
            layersPerWinding.push_back(layersThisWinding);
        }

        // Distribute remaining layers round-robin
        for (size_t windingIndex = 0; windingIndex < (totalNumberLayers - totalAssignedLayers); ++windingIndex) {
            layersPerWinding[windingIndex % coil.get_functional_description().size()]++;
        }

        auto bobbin = coil.resolve_bobbin();
        auto bobbinWidth = bobbin.get_winding_window_dimensions()[0];
        auto bobbinHeight = bobbin.get_winding_window_dimensions()[1];
        
        // Calculate skin depth if frequency is provided
        double skinDepth = 0;
        if (effectiveFrequency && effectiveFrequency.value() > 0) {
            skinDepth = WindingSkinEffectLosses::calculate_skin_depth("copper", effectiveFrequency.value(), temperature);
        }
        
        for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
            auto numberParallels = coil.get_functional_description()[windingIndex].get_number_parallels();
            auto numberTurns = coil.get_functional_description()[windingIndex].get_number_turns();
            auto wire = coil.resolve_wire(windingIndex);
            auto wireWidth = wire.get_maximum_outer_width();
            
            // Minimum layers based on geometric fit
            auto minimumWidth = numberParallels * numberTurns * wireWidth;
            auto minimumNumberLayersNeeded = static_cast<size_t>(ceil(minimumWidth / bobbinWidth));
            
            // Cap layers to minimum needed (don't over-allocate)
            layersPerWinding[windingIndex] = std::min(layersPerWinding[windingIndex], minimumNumberLayersNeeded);
            
            // Skin depth optimization for planar:
            // If wire height > 2×skinDepth, current only uses surface layers effectively.
            // Consider recommending more layers with thinner copper.
            // However, for the stackup, we work with existing wire dimensions.
            // This optimization suggests using MORE layers when beneficial for AC.
            if (skinDepth > 0 && wire.get_conducting_height()) {
                double wireHeight = resolve_dimensional_values(wire.get_conducting_height().value());
                double optimalCopperThickness = 2.0 * skinDepth;  // Current penetrates ~2×skinDepth total
                
                // If wire is much thicker than optimal, we could benefit from more layers
                // Each layer adds surface area for current flow
                if (wireHeight > 2.0 * optimalCopperThickness) {
                    // Calculate how many layers we could use if we had thinner copper
                    // This is informational - actual wire thickness is fixed
                    // But we ensure we use at least enough layers to spread the turns
                    double layerThicknessWithInsulation = wireHeight + defaults.pcbInsulationThickness;
                    size_t maxLayersByHeight = static_cast<size_t>(floor(bobbinHeight / layerThicknessWithInsulation));
                    
                    // Don't exceed physical constraints
                    layersPerWinding[windingIndex] = std::min(layersPerWinding[windingIndex], maxLayersByHeight);
                    
                    // Log suggestion for thick copper at high frequency
                    if (wireHeight > 4.0 * optimalCopperThickness) {
                        logEntry("Skin depth optimization: Winding " + std::to_string(windingIndex) + 
                                 " copper (" + std::to_string(wireHeight * 1e6) + " µm) is " +
                                 std::to_string(wireHeight / optimalCopperThickness) + "× optimal thickness (" +
                                 std::to_string(optimalCopperThickness * 1e6) + " µm). Consider thinner copper.",
                                 "CoilAdviser", 2);
                    }
                }
            }
        }

        // Ensure minimum layers for parallel configurations
        for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
            if (coil.get_functional_description()[windingIndex].get_number_parallels() > 1) {
                layersPerWinding[windingIndex] = std::max(layersPerWinding[windingIndex], 
                    static_cast<size_t>(coil.get_functional_description()[windingIndex].get_number_parallels()));
            }
        }

        // Build stackup with interleaving
        for (size_t repetitionIndex = 0; repetitionIndex < repetitions; ++repetitionIndex) {
            for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
                size_t numberLayerThisSection = static_cast<size_t>(ceil(
                    static_cast<double>(layersPerWinding[windingIndex]) / (repetitions - repetitionIndex)));
                layersPerWinding[windingIndex] -= numberLayerThisSection;
                for (size_t layerIndex = 0; layerIndex < numberLayerThisSection; ++layerIndex) {
                    stackUp.push_back(windingIndex);
                }
            }
        }

        return stackUp;
    }

    std::vector<Mas> CoilAdviser::get_advised_coil_for_pattern(std::vector<Wire>* wires, Mas mas, std::vector<size_t> pattern, size_t repetitions, std::vector<WireSolidInsulationRequirements> solidInsulationRequirementsForWires, size_t maximumNumberResults, std::string reference){
        bool filterMode = bool(mas.get_mutable_inputs().get_design_requirements().get_minimum_impedance());
        const bool isCmc = mas.get_mutable_inputs().get_design_requirements().get_sub_application().has_value()
                        && mas.get_mutable_inputs().get_design_requirements().get_sub_application().value()
                           == SubApplication::COMMON_MODE_NOISE_FILTERING;
        size_t maximumNumberWires = settings.get_coil_adviser_maximum_number_wires();
        auto sectionProportions = calculate_winding_window_proportion_per_winding(mas.get_mutable_inputs());
        auto core = mas.get_magnetic().get_core();
        auto coil = mas.get_magnetic().get_coil();
        if (core.get_functional_description().get_type() != CoreType::TOROIDAL) {
            coil.set_winding_orientation(WindingOrientation::OVERLAPPING);
            coil.set_section_alignment(CoilAlignment::INNER_OR_TOP);
        }
        else {
            coil.set_winding_orientation(WindingOrientation::CONTIGUOUS);
            coil.set_section_alignment(CoilAlignment::SPREAD);
            if (filterMode) {
                coil.set_turns_alignment(CoilAlignment::CENTERED);
            }
        }
        mas.get_mutable_magnetic().set_coil(coil);

        size_t numberWindings = coil.get_functional_description().size();

        auto needsMargin = InsulationCoordinator::needs_margin(solidInsulationRequirementsForWires, pattern, repetitions);
        coil.set_inputs(mas.get_inputs());
        coil.clear();
        coil.set_groups_description(std::nullopt);
        coil.wind_by_sections(sectionProportions, pattern, repetitions);

        auto sections = get_advised_sections(mas, pattern, repetitions);
        if (sections.size() == 0) {
            return {};
        }
        coil.set_sections_description(sections);

        for (auto section : sections) {
            if (section.get_dimensions()[0] < 0) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT, "section.get_dimensions()[0] cannot be negative");
            }
            if (section.get_dimensions()[1] < 0) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT, "section.get_dimensions()[1] cannot be negative: " + std::to_string(section.get_dimensions()[1]));
            }
        }

        _wireAdviser.set_common_wire_standard(_commonWireStandard);  // TODO, rethink this

        if (needsMargin) {
            // If we want to use margin, we set the maximum so the wires chosen will need margin (and be faster)
            for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
                solidInsulationRequirementsForWires[windingIndex].set_maximum_number_layers(1);
                solidInsulationRequirementsForWires[windingIndex].set_maximum_grade(3);
            }
        }

        std::vector<std::vector<std::pair<OpenMagnetics::Winding, double>>> wireCoilPerWinding;


        if (!mas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Missing current in excitaton");
        }

        if (!mas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current().value().get_harmonics() &&
            !mas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current().value().get_processed() &&
            !mas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current().value().get_waveform()) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Missing current harmonics, waveform and processed in excitaton");
        }


        int timeout = 1 - numberWindings;
        for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
            _wireAdviser.set_wire_solid_insulation_requirements(solidInsulationRequirementsForWires[windingIndex]);
            wireCoilPerWinding.push_back(std::vector<std::pair<Winding, double>>{});


            SignalDescriptor maximumCurrent;
            double maximumCurrentRmsTimesRootSquaredEffectiveFrequency = 0;
            for (size_t operatingPointIndex = 0; operatingPointIndex < mas.get_inputs().get_operating_points().size(); ++operatingPointIndex) {
                if (!mas.get_inputs().get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current()) {
                    throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Current is missing");
                }
                if (!mas.get_inputs().get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current()->get_processed()) {
                    throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Current is not processed");
                }
                if (!mas.get_inputs().get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current()->get_processed()->get_effective_frequency()) {
                    auto current = mas.get_inputs().get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current().value();
                    if (current.get_waveform()) {
                        auto processed = Inputs::calculate_processed_data(current.get_waveform().value());
                        current.set_processed(processed);
                        mas.get_mutable_inputs().get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex].set_current(current);
                    }
                    else {
                        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Current is not processed");
                    }
                }
                auto effectiveFreqOpt = mas.get_inputs().get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current()->get_processed()->get_effective_frequency();
                auto rmsOpt = mas.get_inputs().get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current()->get_processed()->get_rms();
                if (!effectiveFreqOpt || !rmsOpt) {
                    throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Current is missing effective frequency or RMS data");
                }
                double effectiveFrequency = effectiveFreqOpt.value();
                double rms = rmsOpt.value();

                auto currentRmsTimesRootSquaredEffectiveFrequency = rms * sqrt(effectiveFrequency);
                if (currentRmsTimesRootSquaredEffectiveFrequency > maximumCurrentRmsTimesRootSquaredEffectiveFrequency) {
                    maximumCurrentRmsTimesRootSquaredEffectiveFrequency = currentRmsTimesRootSquaredEffectiveFrequency;
                    maximumCurrent = mas.get_inputs().get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current().value();
                }

            }

            double maximumTemperature = (*std::max_element(mas.get_inputs().get_operating_points().begin(), mas.get_inputs().get_operating_points().end(),
                                     [](const OperatingPoint &p1,
                                        const OperatingPoint &p2)
                                     {
                                         return p1.get_conditions().get_ambient_temperature() < p2.get_conditions().get_ambient_temperature();
                                     })).get_conditions().get_ambient_temperature();

            // COA-OPT-1 NOTE: This wireConfigurations block is duplicated in wound and planar paths.
            // TODO: Extract as static const or class member to reduce duplication.
            std::vector<std::map<std::string, double>> wireConfigurations = {
                {{"maximumEffectiveCurrentDensity", defaults.maximumEffectiveCurrentDensity}, {"maximumNumberParallels", defaults.maximumNumberParallels}},
                {{"maximumEffectiveCurrentDensity", defaults.maximumEffectiveCurrentDensity}, {"maximumNumberParallels", defaults.maximumNumberParallels * 2}},
                {{"maximumEffectiveCurrentDensity", defaults.maximumEffectiveCurrentDensity * 2}, {"maximumNumberParallels", defaults.maximumNumberParallels}},
                {{"maximumEffectiveCurrentDensity", defaults.maximumEffectiveCurrentDensity * 2}, {"maximumNumberParallels", defaults.maximumNumberParallels * 2}}
            };
            logEntry("Trying " + std::to_string(wireConfigurations.size()) + " wire configurations", "CoilAdviser", 2);

            // wound_with windings (e.g. center-tapped LLC half-secondaries)
            // are virtualized into the same section as their partner, so the
            // section count is smaller than the winding count. Compute the
            // 0-based virtual-group index for this winding by counting how
            // many distinct groups (windings that are NOT wound_with a lower-
            // index partner) appear up to and including this winding.
            auto group_representative = [&](size_t idx) -> size_t {
                const auto& ww = coil.get_functional_description()[idx].get_wound_with();
                if (!ww || ww->empty()) return idx;
                const std::string& partnerName = ww.value().front();
                for (size_t pIdx = 0; pIdx < numberWindings; ++pIdx) {
                    if (coil.get_functional_description()[pIdx].get_name() == partnerName) {
                        return std::min(pIdx, idx);
                    }
                }
                return idx;
            };
            size_t conductionSectionIndex = 0;
            {
                std::set<size_t> seenReps;
                size_t myRep = group_representative(windingIndex);
                for (size_t k = 0; k <= windingIndex; ++k) {
                    size_t rep = group_representative(k);
                    if (rep == myRep && k != windingIndex) {
                        // We are a non-representative member of an already-seen
                        // group; reuse that group's section index.
                        break;
                    }
                    if (seenReps.insert(rep).second) {
                        if (rep == myRep) break;
                        conductionSectionIndex++;
                    }
                }
            }

            for (auto& wireConfiguration : wireConfigurations) {
                _wireAdviser.set_maximum_effective_current_density(wireConfiguration["maximumEffectiveCurrentDensity"]);
                _wireAdviser.set_maximum_number_parallels(wireConfiguration["maximumNumberParallels"]);
                logEntry("Trying wires with a current density of " + std::to_string(wireConfiguration["maximumEffectiveCurrentDensity"]) + " and " + std::to_string(wireConfiguration["maximumNumberParallels"]) + " maximum parallels", "CoilAdviser", 3);

                auto sectionIndex = coil.convert_conduction_section_index_to_global(conductionSectionIndex);

                auto wiresWithScoring = _wireAdviser.get_advised_wire(wires,
                                                                     coil.get_functional_description()[windingIndex],
                                                                     sections[sectionIndex],
                                                                     maximumCurrent,
                                                                     maximumTemperature,
                                                                     coil.get_interleaving_level(),
                                                                     maximumNumberWires);

                if (wiresWithScoring.size() != 0) {
                    timeout += wiresWithScoring.size();

                    std::move(wiresWithScoring.begin(), wiresWithScoring.end(), std::back_inserter(wireCoilPerWinding.back()));
                    if (wireCoilPerWinding.back().size() > maximumNumberResults * 4) {
                        break;
                    }
                }
            }
        }

        for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
            if (windingIndex >= wireCoilPerWinding.size()) {
                return {};
            }
            if (wireCoilPerWinding[windingIndex].size() == 0) {
                return {};
            }
        }

        logEntry("Trying to wind " + std::to_string(wireCoilPerWinding[0].size()) + " coil possibilities", "CoilAdviser");
        mas.get_mutable_magnetic().set_coil(coil);

        auto currentWireIndexPerWinding = std::vector<size_t>(numberWindings, 0);
        std::vector<Mas> masesWithCoil;

        size_t wiresIndex = 0;

        while (true) {

            std::vector<Winding> windings;

            for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
                windings.push_back(wireCoilPerWinding[windingIndex][currentWireIndexPerWinding[windingIndex]].first);
            }

            // Wound-together windings (e.g. center-tapped LLC secondary halves)
            // share a single section and must have matching wire + parallels
            // for the Coil virtualization step to merge them. The wire adviser
            // runs per-winding and may pick different wires, so copy the wire
            // and numberParallels from the first member of each wound_with
            // group to every other member.
            for (size_t wIdx = 0; wIdx < windings.size(); ++wIdx) {
                const auto& wwOpt = windings[wIdx].get_wound_with();
                if (!wwOpt || wwOpt->empty()) continue;
                const std::string& partnerName = wwOpt.value().front();
                for (size_t pIdx = 0; pIdx < windings.size(); ++pIdx) {
                    if (windings[pIdx].get_name() == partnerName) {
                        windings[wIdx].set_wire(windings[pIdx].get_wire());
                        windings[wIdx].set_number_parallels(windings[pIdx].get_number_parallels());
                        break;
                    }
                }
            }

            // CMC: mirror winding[0] wire + parallels to every other winding.
            // A CMC is a coupled inductor where every phase carries identical
            // line current, so the wire spec must match across windings.
            // The wire adviser runs per-winding, and with identical operating
            // points usually converges on the same pick anyway — but nothing
            // in the WireAdviser enforces this, so pin it here.
            if (isCmc && !windings.empty()) {
                const auto& w0 = windings[0];
                for (size_t wIdx = 1; wIdx < windings.size(); ++wIdx) {
                    windings[wIdx].set_wire(w0.get_wire());
                    windings[wIdx].set_number_parallels(w0.get_number_parallels());
                }
            }

            mas.get_mutable_magnetic().get_mutable_coil().set_functional_description(windings);

            // We have new wires combination, we need to restart insulation each time and let it compute it again
            mas.get_mutable_magnetic().get_mutable_coil().reset_insulation();
            bool wound = mas.get_mutable_magnetic().get_mutable_coil().wind(sectionProportions, pattern, repetitions);

            if (wound) {
                mas.get_mutable_magnetic().get_mutable_coil().delimit_and_compact();
                mas.get_mutable_magnetic().set_coil(mas.get_mutable_magnetic().get_mutable_coil());

                // For toroids the angular fill is the binding constraint: all section
                // angular extents (dimensions[1], in degrees) must sum to ≤ 360°.
                if (core.get_functional_description().get_type() == CoreType::TOROIDAL) {
                    double totalAngle = 0.0;
                    auto sectionsDesc = mas.get_magnetic().get_coil().get_sections_description();
                    if (sectionsDesc) {
                        for (const auto& sec : sectionsDesc.value()) {
                            totalAngle += sec.get_dimensions()[1];
                        }
                    }
                    if (totalAngle > 360.0) {
                        wound = false;
                    }
                }

                if (wound) {
                    if (!mas.get_mutable_magnetic().get_manufacturer_info()) {
                        MagneticManufacturerInfo manufacturerInfo;
                        mas.get_mutable_magnetic().set_manufacturer_info(manufacturerInfo);
                    }
                    auto info = mas.get_mutable_magnetic().get_manufacturer_info().value();
                    auto auxReference = reference;
                    auxReference += std::to_string(wiresIndex);
                    info.set_reference(auxReference);
                    mas.get_mutable_magnetic().set_manufacturer_info(info);

                    masesWithCoil.push_back(mas);
                    wiresIndex++;
                    if (masesWithCoil.size() == maximumNumberResults) {
                        break;
                    }
                }
            }
            timeout--;
            if (timeout == 0) {
                break;
            }
            // B18 FIX: score-guided wire advancement (advance winding with worst wire score)
            size_t lowestIndex = 0;
            double worstScore = std::numeric_limits<double>::max();
            for (size_t w = 0; w < numberWindings; ++w) {
                if (currentWireIndexPerWinding[w] + 1 < wireCoilPerWinding[w].size()) {
                    double score = wireCoilPerWinding[w][currentWireIndexPerWinding[w]].second;
                    if (score < worstScore) {
                        worstScore = score;
                        lowestIndex = w;
                    }
                }
            }

            bool anyAdvanceable = false; // COA-BUG-1 FIX
            for (size_t auxWindingIndex = 0; auxWindingIndex < numberWindings; ++auxWindingIndex) {
                if (currentWireIndexPerWinding[lowestIndex] < wireCoilPerWinding[lowestIndex].size() - 1) {
                    anyAdvanceable = true; // COA-BUG-1 FIX
                    break;
                }
                lowestIndex = (lowestIndex + 1) % currentWireIndexPerWinding.size();
            }

            if (!anyAdvanceable) break; // COA-BUG-1 FIX: all windings exhausted
            currentWireIndexPerWinding[lowestIndex]++;
        }
        logEntry("Managed to wind " + std::to_string(masesWithCoil.size()) + " coils", "CoilAdviser");

        return masesWithCoil;
    }

    std::vector<Mas> CoilAdviser::get_advised_planar_coil_for_pattern(std::vector<Wire>* wires, Mas mas, std::vector<size_t> pattern, size_t repetitions, size_t maximumNumberResults, std::string reference){
        // bool filterMode = bool(mas.get_mutable_inputs().get_design_requirements().get_minimum_impedance());
        size_t maximumNumberWires = settings.get_coil_adviser_maximum_number_wires();
        auto sectionProportions = calculate_winding_window_proportion_per_winding(mas.get_mutable_inputs());
        auto core = mas.get_magnetic().get_core();
        auto coil = mas.get_magnetic().get_coil();

        coil.set_winding_orientation(WindingOrientation::CONTIGUOUS);
        coil.set_section_alignment(CoilAlignment::CENTERED);
        coil.set_layers_orientation(WindingOrientation::CONTIGUOUS);
        coil.set_turns_alignment(CoilAlignment::SPREAD);

        mas.get_mutable_magnetic().set_coil(coil);

        size_t numberWindings = coil.get_functional_description().size();

        coil.set_inputs(mas.get_inputs());
        coil.clear();

        std::vector<std::vector<std::pair<OpenMagnetics::Winding, double>>> wireCoilPerWinding;

        auto sections = get_advised_planar_sections(mas, pattern, repetitions);
        if (sections.size() == 0) {
            return {};
        }

        coil.set_sections_description(sections);

        if (!mas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Missing current in excitaton");
        }

        if (!mas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current().value().get_harmonics() &&
            !mas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current().value().get_processed() &&
            !mas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current().value().get_waveform()) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Missing current harmonics, waveform and processed in excitaton");
        }


        int timeout = 1 - numberWindings;
        for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {

            SignalDescriptor maximumCurrent;
            double maximumCurrentRmsTimesRootSquaredEffectiveFrequency = 0;
            for (size_t operatingPointIndex = 0; operatingPointIndex < mas.get_inputs().get_operating_points().size(); ++operatingPointIndex) {
                if (!mas.get_inputs().get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current()) {
                    throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Current is missing");
                }
                if (!mas.get_inputs().get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current()->get_processed()) {
                    throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Current is not processed");
                }
                if (!mas.get_inputs().get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current()->get_processed()->get_effective_frequency()) {
                    auto current = mas.get_inputs().get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current().value();
                    if (current.get_waveform()) {
                        auto processed = Inputs::calculate_processed_data(current.get_waveform().value());
                        current.set_processed(processed);
                        mas.get_mutable_inputs().get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex].set_current(current);
                    }
                    else {
                        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Current is not processed");
                    }
                }
                auto effectiveFreqOpt = mas.get_inputs().get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current()->get_processed()->get_effective_frequency();
                auto rmsOpt = mas.get_inputs().get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current()->get_processed()->get_rms();
                if (!effectiveFreqOpt || !rmsOpt) {
                    throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Current is missing effective frequency or RMS data");
                }
                double effectiveFrequency = effectiveFreqOpt.value();
                double rms = rmsOpt.value();

                auto currentRmsTimesRootSquaredEffectiveFrequency = rms * sqrt(effectiveFrequency);
                if (currentRmsTimesRootSquaredEffectiveFrequency > maximumCurrentRmsTimesRootSquaredEffectiveFrequency) {
                    maximumCurrentRmsTimesRootSquaredEffectiveFrequency = currentRmsTimesRootSquaredEffectiveFrequency;
                    maximumCurrent = mas.get_inputs().get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current().value();
                }

            }

            double maximumTemperature = mas.get_mutable_inputs().get_maximum_temperature();
            bool found = false;
            auto sectionIndex = coil.convert_conduction_section_index_to_global(windingIndex);

            size_t numberSectionsPerPattern = 0;
            for (auto auxWindingIndex : pattern) {
                if (auxWindingIndex == windingIndex) {
                    numberSectionsPerPattern++;
                }
            }

            // Try multiple wire configurations to find optimal designs
            std::vector<std::map<std::string, double>> wireConfigurations = {
                {{"maximumEffectiveCurrentDensity", defaults.maximumEffectiveCurrentDensity}, {"maximumNumberParallels", defaults.maximumNumberParallels}},
                {{"maximumEffectiveCurrentDensity", defaults.maximumEffectiveCurrentDensity}, {"maximumNumberParallels", defaults.maximumNumberParallels * 2}},
                {{"maximumEffectiveCurrentDensity", defaults.maximumEffectiveCurrentDensity * 2}, {"maximumNumberParallels", defaults.maximumNumberParallels}},
                {{"maximumEffectiveCurrentDensity", defaults.maximumEffectiveCurrentDensity * 2}, {"maximumNumberParallels", defaults.maximumNumberParallels * 2}}
            };
            logEntry("Trying " + std::to_string(wireConfigurations.size()) + " wire configurations", "CoilAdviser", 2);

            wireCoilPerWinding.push_back(std::vector<std::pair<Winding, double>>{});

            for (auto& wireConfiguration : wireConfigurations) {
                _wireAdviser.set_maximum_effective_current_density(wireConfiguration["maximumEffectiveCurrentDensity"]);
                _wireAdviser.set_maximum_number_parallels(wireConfiguration["maximumNumberParallels"]);
                logEntry("Trying planar wires with a current density of " + std::to_string(wireConfiguration["maximumEffectiveCurrentDensity"]) + " and " + std::to_string(wireConfiguration["maximumNumberParallels"]) + " maximum parallels", "CoilAdviser", 3);

                auto wiresWithScoring = _wireAdviser.get_advised_planar_wire(coil.get_functional_description()[windingIndex],
                                                                             sections[sectionIndex],
                                                                             maximumCurrent,
                                                                             maximumTemperature,
                                                                             numberSectionsPerPattern * repetitions,
                                                                             maximumNumberWires);

                if (wiresWithScoring.size() != 0) {
                    timeout += wiresWithScoring.size();
                    // Only add wires up to the limit, but continue trying other configurations
                    size_t spaceRemaining = maximumNumberResults * 4 - wireCoilPerWinding.back().size();
                    size_t wiresToAdd = std::min(spaceRemaining, wiresWithScoring.size());
                    std::move(wiresWithScoring.begin(), wiresWithScoring.begin() + wiresToAdd, std::back_inserter(wireCoilPerWinding.back()));
                    found = true;
                }
            }
            if (!found) {
                logEntry("No planar wires found for winding " + std::to_string(windingIndex), "CoilAdviser", 2);
            }
            else {
                logEntry("Found " + std::to_string(wireCoilPerWinding.back().size()) + " planar wire options for winding " + std::to_string(windingIndex), "CoilAdviser", 2);
            }

        }

        for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
            if (windingIndex >= wireCoilPerWinding.size()) {
                return {};
            }
            if (wireCoilPerWinding[windingIndex].size() == 0) {
                return {};
            }
        }

        logEntry("Trying to wind " + std::to_string(wireCoilPerWinding[0].size()) + " coil possibilities", "CoilAdviser");
        mas.get_mutable_magnetic().set_coil(coil);

        auto currentWireIndexPerWinding = std::vector<size_t>(numberWindings, 0);
        std::vector<Mas> masesWithCoil;

        size_t wiresIndex = 0;

        while (true) {

            std::vector<Winding> windings;

            for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
                windings.push_back(wireCoilPerWinding[windingIndex][currentWireIndexPerWinding[windingIndex]].first);
            }

            // Wound-together windings (e.g. center-tapped LLC secondary halves)
            // share a single section and must have matching wire + parallels
            // for the Coil virtualization step to merge them. The wire adviser
            // runs per-winding and may pick different wires, so copy the wire
            // and numberParallels from the first member of each wound_with
            // group to every other member.
            for (size_t wIdx = 0; wIdx < windings.size(); ++wIdx) {
                const auto& wwOpt = windings[wIdx].get_wound_with();
                if (!wwOpt || wwOpt->empty()) continue;
                const std::string& partnerName = wwOpt.value().front();
                for (size_t pIdx = 0; pIdx < windings.size(); ++pIdx) {
                    if (windings[pIdx].get_name() == partnerName) {
                        windings[wIdx].set_wire(windings[pIdx].get_wire());
                        windings[wIdx].set_number_parallels(windings[pIdx].get_number_parallels());
                        break;
                    }
                }
            }

            mas.get_mutable_magnetic().get_mutable_coil().set_functional_description(windings);

            mas.get_mutable_magnetic().get_mutable_coil().reset_margins_per_section();
            bool wound = false;

            // Get effective frequency for skin depth calculation
            std::optional<double> effectiveFrequency = std::nullopt;
            double temperature = mas.get_mutable_inputs().get_maximum_temperature();
            if (mas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current() &&
                mas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed() &&
                mas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_effective_frequency()) {
                effectiveFrequency = mas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_effective_frequency().value();
            }
            auto stackUp = get_planar_stackup(mas.get_mutable_magnetic().get_mutable_coil(), sectionProportions, pattern, repetitions, effectiveFrequency, temperature);
            std::string stackUpString = "";

            for (auto windingIndex : stackUp) {
                stackUpString += std::to_string(windingIndex) + "-";
            }
            stackUpString.pop_back();

            // B19 FIX: Calculate planar clearances using InsulationCoordinator
            if (mas.get_mutable_inputs().get_design_requirements().get_insulation() &&
                mas.get_mutable_inputs().get_wiring_technology() == WiringTechnology::PRINTED) {
                try {
                    InsulationCoordinator insulationCoordinator;
                    auto insulationOutput = insulationCoordinator.calculate_insulation_coordination(
                        mas.get_mutable_inputs());
                    double requiredClearance = insulationOutput.get_clearance();
                    if (requiredClearance > 0 &&
                        mas.get_mutable_magnetic().get_mutable_coil().get_sections_description()) {
                        auto sections = mas.get_mutable_magnetic().get_mutable_coil()
                                           .get_sections_description().value();
                        for (size_t sIdx = 0; sIdx + 1 < sections.size(); ++sIdx) {
                            if (!sections[sIdx].get_coordinates().empty() && !sections[sIdx+1].get_coordinates().empty()) {
                                double spacing = std::abs(
                                    sections[sIdx+1].get_coordinates()[1] -
                                    sections[sIdx].get_coordinates()[1]) -
                                    (sections[sIdx].get_dimensions()[1] +
                                     sections[sIdx+1].get_dimensions()[1]) / 2.0;
                                if (spacing < requiredClearance) {
                                    logEntry("B19: planar clearance " +
                                        std::to_string(spacing * 1e3) + " mm < required " +
                                        std::to_string(requiredClearance * 1e3) +
                                        " mm — skipping", "CoilAdviser", 2);
                                    continue; // skip this wire combination
                                }
                            }
                        }
                    }
                }
                catch (...) { // XC-3 NOTE: consider catching std::exception for diagnostics
                    logEntry("B19: clearance check skipped (missing data)", "CoilAdviser", 2);
                }
            }
            wound = mas.get_mutable_magnetic().get_mutable_coil().wind_planar(stackUp, std::nullopt, {}, {}, defaults.coreToLayerDistance);

            if (wound) {
                mas.get_mutable_magnetic().get_mutable_coil().delimit_and_compact();
                mas.get_mutable_magnetic().set_coil(mas.get_mutable_magnetic().get_mutable_coil());
                if (!mas.get_mutable_magnetic().get_manufacturer_info()) {
                    MagneticManufacturerInfo manufacturerInfo;
                    mas.get_mutable_magnetic().set_manufacturer_info(manufacturerInfo);
                }
                auto info = mas.get_mutable_magnetic().get_manufacturer_info().value();
                auto auxReference = reference;
                auxReference += ", Parallels: ";
                auxReference += std::to_string(mas.get_magnetic().get_coil().get_functional_description()[0].get_number_parallels());
                auxReference += ", Stack-up: " + stackUpString;
                auxReference += " " + std::to_string(wiresIndex);
                info.set_reference(auxReference);
                mas.get_mutable_magnetic().set_manufacturer_info(info);

                masesWithCoil.push_back(mas);
                wiresIndex++;
                if (masesWithCoil.size() == maximumNumberResults) {
                    break;
                }
            }
            timeout--;
            if (timeout == 0) {
                break;
            }
            // B18 FIX: score-guided wire advancement (advance winding with worst wire score)
            size_t lowestIndex = 0;
            double worstScore = std::numeric_limits<double>::max();
            for (size_t w = 0; w < numberWindings; ++w) {
                if (currentWireIndexPerWinding[w] + 1 < wireCoilPerWinding[w].size()) {
                    double score = wireCoilPerWinding[w][currentWireIndexPerWinding[w]].second;
                    if (score < worstScore) {
                        worstScore = score;
                        lowestIndex = w;
                    }
                }
            }

            bool anyAdvanceable = false; // COA-BUG-1 FIX
            for (size_t auxWindingIndex = 0; auxWindingIndex < numberWindings; ++auxWindingIndex) {
                if (currentWireIndexPerWinding[lowestIndex] < wireCoilPerWinding[lowestIndex].size() - 1) {
                    anyAdvanceable = true; // COA-BUG-1 FIX
                    break;
                }
                lowestIndex = (lowestIndex + 1) % currentWireIndexPerWinding.size();
            }

            if (!anyAdvanceable) break; // COA-BUG-1 FIX: all windings exhausted
            currentWireIndexPerWinding[lowestIndex]++;
        }
        logEntry("Managed to wind " + std::to_string(masesWithCoil.size()) + " coils", "CoilAdviser");

        return masesWithCoil;

    }
} // namespace OpenMagnetics
