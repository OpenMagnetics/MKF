#include "processors/Inputs.h"
#include <source_location>
#include "advisers/MagneticAdviser.h"
#include "physical_models/Impedance.h"
#include "processors/MagneticSimulator.h"
#include "support/Painter.h"
#include <magic_enum_utility.hpp>
#include "support/Logger.h"


namespace OpenMagnetics {

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

std::vector<std::pair<Mas, double>> MagneticAdviser::get_advised_magnetic(Inputs inputs, size_t maximumNumberResults) {
    return get_advised_magnetic(inputs, _defaultCustomMagneticFilterFlow, maximumNumberResults);
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
    
    // Store original toroid setting for potential retry
    bool toroidsOriginallyEnabled = settings.get_use_toroidal_cores();

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
        logEntry("High voltage requirements (" + std::to_string(int(maxVoltage)) + "V) detected. Disabling toroidal cores for better results.", "MagneticAdviser", 1);
        settings.set_use_toroidal_cores(false);
        toroidsOriginallyEnabled = false; // Don't retry since we proactively disabled
    }

    if (get_application() == Application::INTERFERENCE_SUPPRESSION) {
        settings.set_use_toroidal_cores(true);
        settings.set_use_only_cores_in_stock(false);
        settings.set_use_concentric_cores(false);
    }

    if (coreDatabase.empty()) {
        load_cores();
    }
    if (wireDatabase.empty()) {
        load_wires();
    }

    bool previousCoilIncludeAdditionalCoordinates = settings.get_coil_include_additional_coordinates();
    settings.set_coil_include_additional_coordinates(false);

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
    while (coresWound < expectedWoundCores && whileIteration < maxWhileIterations && evaluatedCores.size() < maxEvaluatedCores) {
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
                mas = magneticSimulator.simulate(mas);

                processedCoils++;

                masData.push_back(mas);
                if (processedCoils >= size_t(ceil(maximumNumberResults * 0.5))) {
                    usedNumberSectionsAndMargin.push_back(numberSectionsAndMarginCombination);
                    break;
                }
            }
            if (coresWound >= expectedWoundCores) {
                break;
            }
        }
    }

    logEntry("Found " + std::to_string(masData.size()) + " magnetics", "MagneticAdviser", 2);
    auto masMagneticsWithScoring = score_magnetics(masData, filterFlow);

    sort(masMagneticsWithScoring.begin(), masMagneticsWithScoring.end(), [](std::pair<Mas, double>& b1, std::pair<Mas, double>& b2) {
        return b1.second > b2.second;
    });

    if (masMagneticsWithScoring.size() > maximumNumberResults) {
        masMagneticsWithScoring = std::vector<std::pair<Mas, double>>(masMagneticsWithScoring.begin(), masMagneticsWithScoring.end() - (masMagneticsWithScoring.size() - maximumNumberResults));
    }

    // Retry without toroids if toroids were enabled but no results found
    if (masMagneticsWithScoring.empty() && toroidsOriginallyEnabled) {
        logEntry("No magnetics found with toroids enabled. Retrying without toroids...", "MagneticAdviser", 1);
        settings.set_use_toroidal_cores(false);
        // Reset evaluated cores to allow re-evaluation
        evaluatedCores.clear();
        coresWound = 0;
        masData.clear();
        
        // Retry the core search
        whileIteration = 0;
        requestedCores = expectedWoundCores;
        previouslyObtainedCores = SIZE_MAX;
        while (coresWound < expectedWoundCores && whileIteration < maxWhileIterations && evaluatedCores.size() < maxEvaluatedCores) {
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
                        if (processedCoils >= size_t(ceil(maximumNumberResults * 0.5))) {
                            break;
                        }
                    }
                }
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

    settings.set_coil_include_additional_coordinates(previousCoilIncludeAdditionalCoordinates);
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

std::vector<std::pair<Mas, double>> MagneticAdviser::get_advised_magnetic(Inputs inputs, std::map<std::string, Magnetic> catalogueMagnetics, std::vector<MagneticFilterOperation> filterFlow, size_t maximumNumberResults, bool strict) {
    std::vector<Mas> catalogueMagneticsWithInputs;
    for (auto [reference, magnetic] : catalogueMagnetics) {
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

std::vector<std::pair<Mas, double>> MagneticAdviser::get_advised_magnetic(std::vector<Mas> catalogueMagneticsWithInputs, std::vector<MagneticFilterOperation> filterFlow, size_t maximumNumberResults, bool strict) {

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
            }
            catch (...) {
                validMagnetic = false;
                break;
            }
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
        std::vector<Outputs> outputs;
        auto inputs = mas.get_inputs();
        auto magnetic = mas.get_magnetic();
        bool valid = true;
        for (auto filterConfiguration : filterFlow) {
            MagneticFilters filterEnum = filterConfiguration.get_filter();
        
            try {
                auto [valid, scoring] = _filters[filterEnum]->evaluate_magnetic(&magnetic, &inputs, &outputs);
                add_scoring(magnetic.get_reference(), filterEnum, scoring);
            }
            catch (...) {
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
                mas = magneticSimulator.simulate(mas, true);
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
            
            if (coreLosses.get_temperature()) {
                text += "\tCore temperature: " + std::to_string(coreLosses.get_temperature().value()) + "\n";
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
    std::cout << text << std::endl;
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