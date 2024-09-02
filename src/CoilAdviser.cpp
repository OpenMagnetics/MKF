#include "InputsWrapper.h"
#include "WireAdviser.h"
#include "CoilAdviser.h"
#include "Defaults.h"
#include "Models.h"
#include "Insulation.h"
#include "Painter.h"
#include "Settings.h"


namespace OpenMagnetics {
    bool needs_margin(std::vector<WireSolidInsulationRequirements> combinationSolidInsulationRequirementsForWires, std::vector<size_t> pattern, size_t repetitions) {
        for(size_t index = 0; index < pattern.size(); ++index) {
            size_t leftIndex = pattern[index];
            size_t rightIndex;
            if (leftIndex == pattern.size() - 1) {
                if (repetitions == 0) {
                    break;
                }
                else {
                    rightIndex = pattern[0];
                }
            }
            else {
                rightIndex = pattern[leftIndex + 1];
            }
            int64_t numberLayers = 0;
            int64_t numberGrades = 0;
            if (combinationSolidInsulationRequirementsForWires[leftIndex].get_minimum_grade()) {
                numberGrades += combinationSolidInsulationRequirementsForWires[leftIndex].get_minimum_grade().value();
            }
            if (combinationSolidInsulationRequirementsForWires[leftIndex].get_minimum_number_layers()) {
                numberLayers += combinationSolidInsulationRequirementsForWires[leftIndex].get_minimum_number_layers().value();
            }
            if (combinationSolidInsulationRequirementsForWires[rightIndex].get_minimum_grade()) {
                numberGrades += combinationSolidInsulationRequirementsForWires[rightIndex].get_minimum_grade().value();
            }
            if (combinationSolidInsulationRequirementsForWires[rightIndex].get_minimum_number_layers()) {
                numberLayers += combinationSolidInsulationRequirementsForWires[rightIndex].get_minimum_number_layers().value();
            }

            if (numberLayers < 3 && numberGrades < 4) {
                return true;
            }
        }
        return false;
    }

    std::vector<double> calculate_winding_window_proportion_per_winding(InputsWrapper& inputs) {
        auto averagePowerPerWinding = std::vector<double>(inputs.get_operating_points()[0].get_excitations_per_winding().size(), 0);

        for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
            for (size_t windingIndex = 0; windingIndex < inputs.get_operating_points()[operatingPointIndex].get_excitations_per_winding().size(); ++windingIndex) {

                if (!inputs.get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_voltage()) {
                    auto excitation = inputs.get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex];
                    inputs.get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex].set_voltage(InputsWrapper::calculate_induced_voltage(excitation, resolve_dimensional_values(inputs.get_design_requirements().get_magnetizing_inductance())));
                }

                auto power = InputsWrapper::calculate_instantaneous_power(inputs.get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex]);
                averagePowerPerWinding[windingIndex] += power;
            }
        }

        double sumValue = std::reduce(averagePowerPerWinding.begin(), averagePowerPerWinding.end());
        for (size_t windingIndex = 0; windingIndex < averagePowerPerWinding.size(); ++windingIndex) {
            averagePowerPerWinding[windingIndex] = std::max(0.05, averagePowerPerWinding[windingIndex] / sumValue);
        }

        std::vector<double> proportions;
        sumValue = std::reduce(averagePowerPerWinding.begin(), averagePowerPerWinding.end());
        for (size_t windingIndex = 0; windingIndex < averagePowerPerWinding.size(); ++windingIndex) {
            proportions.push_back(averagePowerPerWinding[windingIndex] / sumValue);
        }

        return proportions;
    }


    std::vector<IsolationSide> get_isolation_sides(InputsWrapper& inputs) {
        if (!inputs.get_design_requirements().get_isolation_sides()) {
            std::vector<IsolationSide> isolationSides = {IsolationSide::PRIMARY};
            for (size_t windingIndex = 1; windingIndex < inputs.get_mutable_design_requirements().get_turns_ratios().size() + 1; ++windingIndex) {
                // isolationSides.push_back(IsolationSide::SECONDARY);
                isolationSides.push_back(get_isolation_side_from_index(windingIndex));
            }
            inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        }

        std::vector<IsolationSide> isolationSidesUsed;
        std::vector<IsolationSide> isolationSidesFromRequirements = inputs.get_design_requirements().get_isolation_sides().value();
        for (auto isolationSide : isolationSidesFromRequirements) {
            if (std::find(isolationSidesUsed.begin(), isolationSidesUsed.end(), isolationSide) == isolationSidesUsed.end()) {
                isolationSidesUsed.push_back(isolationSide);
            }
        }

        return isolationSidesUsed;
    }

    std::vector<std::vector<size_t>> get_patterns(InputsWrapper& inputs, CoreType coreType) {
        auto isolationSidesRequired = get_isolation_sides(inputs);

        if (!inputs.get_design_requirements().get_isolation_sides()) {
            throw std::runtime_error("Missing isolation sides requirement");
        }

        auto isolationSidesRequirement = inputs.get_design_requirements().get_isolation_sides().value();

        std::vector<std::vector<size_t>> sectionPatterns;
        for(size_t i = 0; i < tgamma(isolationSidesRequired.size() + 1) / 2; ++i) {
            std::vector<size_t> sectionPattern;
            for (auto isolationSide : isolationSidesRequired) {
                for (size_t windingIndex = 0; windingIndex < inputs.get_mutable_design_requirements().get_turns_ratios().size() + 1; ++windingIndex) {
                    if (isolationSidesRequirement[windingIndex] == isolationSide) {
                        sectionPattern.push_back(windingIndex);
                    }
                }
            }
            sectionPatterns.push_back(sectionPattern);
            if (sectionPatterns.size() > Defaults().maximumCoilPattern) {
                break;
            }

            std::next_permutation(isolationSidesRequired.begin(), isolationSidesRequired.end());
        }

        if (coreType == CoreType::TOROIDAL) {
            // We remove the last combination as in toroids they go around
            size_t elementsToKeep = std::max(size_t(1), isolationSidesRequired.size() - 1);
            sectionPatterns = std::vector<std::vector<size_t>>(sectionPatterns.begin(), sectionPatterns.end() - (sectionPatterns.size() - elementsToKeep));
        }

        return sectionPatterns;
    }


    std::vector<size_t> get_repetitions(InputsWrapper& inputs, CoreType coreType) {
        if (inputs.get_design_requirements().get_turns_ratios().size() == 0 || coreType == CoreType::TOROIDAL) {
            return {1};  // hardcoded
        }
        if (inputs.get_design_requirements().get_leakage_inductance()) {
            return {2, 1};  // hardcoded        
        }
        else{
            return {1, 2};  // hardcoded        
        }
    }


    std::vector<MasWrapper> CoilAdviser::get_advised_coil(MasWrapper mas, size_t maximumNumberResults){
        auto settings = Settings::GetInstance();
        if (wireDatabase.empty()) {
            load_wires();
        }
        std::string jsonLine;
        std::vector<WireWrapper> wires;
        for (const auto& [key, wire] : wireDatabase) {
            if ((settings->get_wire_adviser_include_planar() || wire.get_type() != WireType::PLANAR) &&
                (settings->get_wire_adviser_include_foil() || wire.get_type() != WireType::FOIL) &&
                (settings->get_wire_adviser_include_rectangular() || wire.get_type() != WireType::RECTANGULAR) &&
                (settings->get_wire_adviser_include_litz() || wire.get_type() != WireType::LITZ) &&
                (settings->get_wire_adviser_include_round() || wire.get_type() != WireType::ROUND)) {
                wires.push_back(wire);
            }
        }
        return get_advised_coil(&wires, mas, maximumNumberResults);
    }

    std::pair<std::vector<size_t>, size_t> check_integrity(std::vector<size_t> pattern, size_t repetitions, CoilWrapper coil) {
        bool needsMerge = false;
        for (auto winding : coil.get_functional_description()) {
            // TODO expand for more than one winding per layer
            auto numberPhysicalTurns = winding.get_number_turns() * winding.get_number_parallels();
            if (numberPhysicalTurns < repetitions) {
                needsMerge = true;
            }
        }

        std::vector<size_t> newPattern;
        if (needsMerge) {
            for (size_t repetition = 1; repetition <= repetitions; ++repetition) {
                for (auto windingIndex : pattern) {
                    auto winding = coil.get_functional_description()[windingIndex];
                    auto numberPhysicalTurns = winding.get_number_turns() * winding.get_number_parallels();
                    if (numberPhysicalTurns >= repetition) {
                        newPattern.push_back(windingIndex);
                    }
                }
            }
            return {newPattern, 1};
        }
        return {pattern, repetitions};
    }

    std::vector<MasWrapper> CoilAdviser::get_advised_coil(std::vector<WireWrapper>* wires, MasWrapper mas, size_t maximumNumberResults){
        auto core = mas.get_magnetic().get_core();
        auto coreType = core.get_functional_description().get_type();

        auto inputs = mas.get_mutable_inputs();
        auto patterns = get_patterns(inputs, coreType);
        auto repetitions = get_repetitions(inputs, coreType);
        mas.set_inputs(inputs);

        size_t maximumNumberResultsPerPattern = std::max(2.0, ceil(maximumNumberResults / (patterns.size() * repetitions.size())));

        std::vector<MasWrapper> masMagneticsWithCoil;
        for (auto repetition : repetitions) {
            for (auto pattern : patterns) {
                auto aux = check_integrity(pattern, repetition, mas.get_magnetic().get_coil());
                pattern = aux.first;
                repetition = aux.second;
                auto combinationsSolidInsulationRequirementsForWires = get_solid_insulation_requirements_for_wires(mas.get_mutable_inputs(), pattern, repetition);
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
                    if (needs_margin(solidInsulationRequirementsForWires, pattern, repetition)) {
                        reference += ", Margin Taped";
                    }
                    else {
                        reference += ", Wire Insulated";
                    }
                    reference += " " + std::to_string(insulationIndex);



                    auto resultsPerPattern = get_advised_coil_for_pattern(wires, mas, pattern, repetition, solidInsulationRequirementsForWires, maximumNumberResultsPerPattern, reference);

                    std::move(resultsPerPattern.begin(), resultsPerPattern.end(), std::back_inserter(masMagneticsWithCoil));
                }
            }
        }

        return masMagneticsWithCoil;
    }

    WireSolidInsulationRequirements get_requirements_for_functional() {
        WireSolidInsulationRequirements wireSolidInsulationRequirements;
        wireSolidInsulationRequirements.set_minimum_grade(1);
        wireSolidInsulationRequirements.set_minimum_number_layers(1);
        wireSolidInsulationRequirements.set_minimum_breakdown_voltage(0);
        return wireSolidInsulationRequirements;
    }

    WireSolidInsulationRequirements get_requirements_for_basic(double withstandVoltage, bool canFullyInsulatedWireBeUsed) {
        WireSolidInsulationRequirements wireSolidInsulationRequirements;
        if (canFullyInsulatedWireBeUsed) {
            wireSolidInsulationRequirements.set_minimum_grade(3);
        }
        wireSolidInsulationRequirements.set_minimum_number_layers(1);
        wireSolidInsulationRequirements.set_minimum_breakdown_voltage(withstandVoltage);
        return wireSolidInsulationRequirements;
    }

    WireSolidInsulationRequirements get_requirements_for_reinforced(double withstandVoltage, bool canFullyInsulatedWireBeUsed) {
        WireSolidInsulationRequirements wireSolidInsulationRequirements;
        if (canFullyInsulatedWireBeUsed) {
            wireSolidInsulationRequirements.set_minimum_grade(3);
        }
        wireSolidInsulationRequirements.set_minimum_number_layers(3);
        wireSolidInsulationRequirements.set_minimum_breakdown_voltage(withstandVoltage);
        return wireSolidInsulationRequirements;
    }

    std::vector<std::vector<WireSolidInsulationRequirements>> remove_combination_that_need_margin(std::vector<std::vector<WireSolidInsulationRequirements>> combinationsSolidInsulationRequirementsForWires, std::vector<size_t> pattern, size_t repetitions) {
        std::vector<std::vector<WireSolidInsulationRequirements>> combinationsSolidInsulationRequirementsForWiresWithoutMargin;
        for (auto combinationSolidInsulationRequirementsForWires : combinationsSolidInsulationRequirementsForWires) {
            bool needsMargin = needs_margin(combinationSolidInsulationRequirementsForWires, pattern, repetitions);
            if (!needsMargin) {
                combinationsSolidInsulationRequirementsForWiresWithoutMargin.push_back(combinationSolidInsulationRequirementsForWires);
            }
        }
        return combinationsSolidInsulationRequirementsForWiresWithoutMargin;
    }

    std::vector<std::vector<WireSolidInsulationRequirements>> CoilAdviser::get_solid_insulation_requirements_for_wires(InputsWrapper& inputs, std::vector<size_t> pattern, size_t repetitions) {

        auto isolationSidesRequired = get_isolation_sides(inputs);
        auto withstandVoltage = InsulationCoordinator().calculate_withstand_voltage(inputs);
        size_t numberWindings = inputs.get_design_requirements().get_turns_ratios().size() + 1;
        std::vector<std::vector<WireSolidInsulationRequirements>> combinationsSolidInsulationRequirementsForWires;

        if (!inputs.get_design_requirements().get_insulation()) {
            std::vector<WireSolidInsulationRequirements> solidInsulationForWires;
            for(size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
                solidInsulationForWires.push_back(get_requirements_for_functional());
            }
            combinationsSolidInsulationRequirementsForWires.push_back(solidInsulationForWires);
        }
        else {

            auto insulationType = inputs.get_insulation_type();
            bool canFullyInsulatedWireBeUsed = InsulationCoordinator::can_fully_insulated_wire_be_used(inputs);
            auto isolationSidePerWinding = inputs.get_design_requirements().get_isolation_sides().value();
            auto settings = Settings::GetInstance();
            bool allowMarginTape = settings->get_coil_allow_margin_tape();

            // If we don't want margin tape we have to increase to insulation type to DOUBLE:
            if ((insulationType == InsulationType::BASIC || insulationType == InsulationType::SUPPLEMENTARY)) {
                if (!allowMarginTape) {
                    insulationType = InsulationType::DOUBLE;
                }
            }

            // Unless we cannot use margin tape, we always add the option of the wires not complying with anything
            if (allowMarginTape) {
                std::vector<WireSolidInsulationRequirements> solidInsulationForWires;
                for(size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
                    solidInsulationForWires.push_back(get_requirements_for_functional());
                }
                combinationsSolidInsulationRequirementsForWires.push_back(solidInsulationForWires);
            }

            if (insulationType != InsulationType::FUNCTIONAL) {

                if ((insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE)) {
                    // Reinforced or single winding inductors must have all the isolation needed in one later or coating
                    for (auto isolationSideToRemoveInsulation : isolationSidesRequired) {
                        std::vector<WireSolidInsulationRequirements> solidInsulationForWires;
                        for (auto isolationSidePerWinding : isolationSidePerWinding) {
                            if (isolationSidePerWinding == isolationSideToRemoveInsulation) {
                                solidInsulationForWires.push_back(get_requirements_for_functional());
                            }
                            else {
                                solidInsulationForWires.push_back(get_requirements_for_reinforced(withstandVoltage, canFullyInsulatedWireBeUsed));
                            }
                        }
                        combinationsSolidInsulationRequirementsForWires.push_back(solidInsulationForWires);
                    }
                }

                if (insulationType != InsulationType::REINFORCED) {
                    // If insulation is Double, we need to get the withstand voltage for basic, not double
                    if (insulationType == InsulationType::DOUBLE) {
                        auto insulation = inputs.get_mutable_design_requirements().get_insulation().value();
                        insulation.set_insulation_type(InsulationType::BASIC);
                        inputs.get_mutable_design_requirements().set_insulation(insulation);
                        withstandVoltage = InsulationCoordinator().calculate_withstand_voltage(inputs);
                    }

                    if (insulationType == InsulationType::DOUBLE || isolationSidesRequired.size() == 1) {
                        // Additionally, Double can be composed of several steps or parts, as long as each reachs Basic and Supplementary insulation.

                        std::vector<WireSolidInsulationRequirements> solidInsulationForWires;
                        for (size_t i = 0; i < isolationSidePerWinding.size(); ++i) {
                            solidInsulationForWires.push_back(get_requirements_for_basic(withstandVoltage, canFullyInsulatedWireBeUsed));
                        }
                        combinationsSolidInsulationRequirementsForWires.push_back(solidInsulationForWires);
                    }

                    if (isolationSidesRequired.size() > 1) {
                        for (auto isolationSideToRemoveInsulation : isolationSidesRequired) {
                            std::vector<WireSolidInsulationRequirements> solidInsulationForWires;
                            for (auto isolationSidePerWinding : isolationSidePerWinding) {
                                if (isolationSidePerWinding == isolationSideToRemoveInsulation) {
                                    solidInsulationForWires.push_back(get_requirements_for_functional());
                                }
                                else {
                                    solidInsulationForWires.push_back(get_requirements_for_basic(withstandVoltage, canFullyInsulatedWireBeUsed));
                                }
                            }
                            combinationsSolidInsulationRequirementsForWires.push_back(solidInsulationForWires);
                        }
                    }
                }
            } 

            if (!allowMarginTape) {
                combinationsSolidInsulationRequirementsForWires = remove_combination_that_need_margin(combinationsSolidInsulationRequirementsForWires, pattern, repetitions);
            }
        }

        return combinationsSolidInsulationRequirementsForWires;
    }

    std::vector<Section> CoilAdviser::get_advised_sections(MasWrapper mas, std::vector<size_t> pattern, size_t repetitions){
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
        }

        coil.set_strict(false);
        coil.set_inputs(mas.get_inputs());
        coil.calculate_insulation(true);
        coil.wind_by_sections(sectionProportions, pattern, repetitions);
        coil.delimit_and_compact();
        coil.set_strict(true);

        auto sections = coil.get_sections_description().value();
        return sections;
    }

    std::vector<MasWrapper> CoilAdviser::get_advised_coil_for_pattern(std::vector<WireWrapper>* wires, MasWrapper mas, std::vector<size_t> pattern, size_t repetitions, std::vector<WireSolidInsulationRequirements> solidInsulationRequirementsForWires, size_t maximumNumberResults, std::string reference){
        auto settings = Settings::GetInstance();
        size_t maximumNumberWires = settings->get_coil_adviser_maximum_number_wires();
        auto defaults = Defaults();
        auto sectionProportions = calculate_winding_window_proportion_per_winding(mas.get_mutable_inputs());
        auto core = mas.get_magnetic().get_core();
        auto coil = mas.get_magnetic().get_coil();
        size_t numberWindings = coil.get_functional_description().size();

        auto needsMargin = needs_margin(solidInsulationRequirementsForWires, pattern, repetitions);
        coil.set_inputs(mas.get_inputs());
        coil.wind_by_sections(sectionProportions, pattern, repetitions);

        auto sections = get_advised_sections(mas, pattern, repetitions);
        coil.set_sections_description(sections);

        for (auto section : sections) {
            if (section.get_dimensions()[0] < 0) {
                throw std::runtime_error("section.get_dimensions()[0] cannot be negative");
            }
            if (section.get_dimensions()[1] < 0) {
                throw std::runtime_error("section.get_dimensions()[1] cannot be negative: " + std::to_string(section.get_dimensions()[1]));
            }
        }

        if (needsMargin) {
            // If we want to use margin, we set the maximum so the wires chosen will need margin (and be faster)
            for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
                solidInsulationRequirementsForWires[windingIndex].set_maximum_number_layers(1);
                solidInsulationRequirementsForWires[windingIndex].set_maximum_grade(3);
            }
        }

        std::vector<std::vector<std::pair<CoilFunctionalDescription, double>>> wireCoilPerWinding;


        if (!mas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()) {
            throw std::runtime_error("Missing current in excitaton");
        }

        if (!mas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current().value().get_harmonics() &&
            !mas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current().value().get_processed() &&
            !mas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current().value().get_waveform()) {
            throw std::runtime_error("Missing current harmonics, waveform and processed in excitaton");
        }


        int timeout = 1 - numberWindings;
        for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
            _wireAdviser.set_wire_solid_insulation_requirements(solidInsulationRequirementsForWires[windingIndex]);

            SignalDescriptor maximumCurrent;
            double maximumCurrentRmsTimesRootSquaredEffectiveFrequency = 0;
            for (size_t operatingPointIndex = 0; operatingPointIndex < mas.get_inputs().get_operating_points().size(); ++operatingPointIndex) {
                if (!mas.get_inputs().get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current()->get_processed()) {
                    throw std::runtime_error("Current is not processed");
                }
                double effectiveFrequency = mas.get_inputs().get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current()->get_processed()->get_effective_frequency().value();
                double rms = mas.get_inputs().get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current()->get_processed()->get_rms().value();

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
            std::vector<std::map<std::string, double>> wireConfigurations = {
                {{"maximumEffectiveCurrentDensity", defaults.maximumEffectiveCurrentDensity}, {"maximumNumberParallels", defaults.maximumNumberParallels}},
                {{"maximumEffectiveCurrentDensity", defaults.maximumEffectiveCurrentDensity}, {"maximumNumberParallels", defaults.maximumNumberParallels * 2}},
                {{"maximumEffectiveCurrentDensity", defaults.maximumEffectiveCurrentDensity * 2}, {"maximumNumberParallels", defaults.maximumNumberParallels}},
                {{"maximumEffectiveCurrentDensity", defaults.maximumEffectiveCurrentDensity * 2}, {"maximumNumberParallels", defaults.maximumNumberParallels * 2}}
            };
            bool found = false;
            for (auto& wireConfiguration : wireConfigurations) {
                _wireAdviser.set_maximum_effective_current_density(wireConfiguration["maximumEffectiveCurrentDensity"]);
                _wireAdviser.set_maximum_number_parallels(wireConfiguration["maximumNumberParallels"]);

                auto sectionIndex = coil.convert_conduction_section_index_to_global(windingIndex);

                auto wiresWithScoring = _wireAdviser.get_advised_wire(wires,
                                                                     coil.get_functional_description()[windingIndex],
                                                                     sections[sectionIndex],
                                                                     maximumCurrent,
                                                                     maximumTemperature,
                                                                     coil.get_interleaving_level(),
                                                                     maximumNumberWires);

                if (wiresWithScoring.size() != 0) {
                    timeout += wiresWithScoring.size();

                    wireCoilPerWinding.push_back(wiresWithScoring);
                    found = true;
                    break;
                }
            }
            if (!found) {
                wireCoilPerWinding.push_back(std::vector<std::pair<CoilFunctionalDescription, double>>{});
            }
        }
        mas.get_mutable_magnetic().set_coil(coil);

        for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
            if (windingIndex >= wireCoilPerWinding.size()) {
                return {};
            }
            if (wireCoilPerWinding[windingIndex].size() == 0) {
                return {};
            }
        }

        auto currentWireIndexPerWinding = std::vector<size_t>(numberWindings, 0);
        std::vector<MasWrapper> masMagneticsWithCoil;

        size_t wiresIndex = 0;

        while (true) {

            std::vector<CoilFunctionalDescription> coilFunctionalDescription;

            for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
                coilFunctionalDescription.push_back(wireCoilPerWinding[windingIndex][currentWireIndexPerWinding[windingIndex]].first);
            }
            mas.get_mutable_magnetic().get_mutable_coil().set_functional_description(coilFunctionalDescription);

            mas.get_mutable_magnetic().get_mutable_coil().reset_margins_per_section();

            bool wound = mas.get_mutable_magnetic().get_mutable_coil().wind(sectionProportions, pattern, repetitions);

            if (wound) {
                mas.get_mutable_magnetic().get_mutable_coil().delimit_and_compact();
                mas.get_mutable_magnetic().set_coil(mas.get_mutable_magnetic().get_mutable_coil());
                if (!mas.get_mutable_magnetic().get_manufacturer_info()) {
                    MagneticManufacturerInfo manufacturerInfo;
                    mas.get_mutable_magnetic().set_manufacturer_info(manufacturerInfo);
                }
                auto info = mas.get_mutable_magnetic().get_manufacturer_info().value();
                auto auxReference = reference;
                auxReference += std::to_string(wiresIndex);
                info.set_reference(auxReference);
                mas.get_mutable_magnetic().set_manufacturer_info(info);

                masMagneticsWithCoil.push_back(mas);
                wiresIndex++;
                if (masMagneticsWithCoil.size() == maximumNumberResults) {
                    break;
                }
            }
            timeout--;
            if (timeout == 0) {
                break;
            }
            auto lowestIndex = std::distance(std::begin(currentWireIndexPerWinding), std::min_element(std::begin(currentWireIndexPerWinding), std::end(currentWireIndexPerWinding)));

            for (size_t auxWindingIndex = 0; auxWindingIndex < numberWindings; ++auxWindingIndex) {
                if (currentWireIndexPerWinding[lowestIndex] < wireCoilPerWinding[lowestIndex].size() - 1) {
                    break;
                }
                lowestIndex = (lowestIndex + 1) % currentWireIndexPerWinding.size();
            }

            currentWireIndexPerWinding[lowestIndex]++;
        }

        return masMagneticsWithCoil;

    }
} // namespace OpenMagnetics
