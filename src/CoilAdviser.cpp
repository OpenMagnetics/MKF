#include "InputsWrapper.h"
#include "WireAdviser.h"
#include "CoilAdviser.h"
#include "Defaults.h"
#include "Models.h"
#include "Insulation.h"
#include "Painter.h"


namespace OpenMagnetics {

    std::vector<double> calculate_winding_window_proportion_per_winding(InputsWrapper inputs) {
        auto averagePowerPerWinding = std::vector<double>(inputs.get_operating_points()[0].get_excitations_per_winding().size(), 0);

        for (size_t operationPointIndex = 0; operationPointIndex < inputs.get_operating_points().size(); ++operationPointIndex) {
            for (size_t windingIndex = 0; windingIndex < inputs.get_operating_points()[operationPointIndex].get_excitations_per_winding().size(); ++windingIndex) {
                auto power = InputsWrapper::calculate_instantaneous_power(inputs.get_operating_points()[operationPointIndex].get_excitations_per_winding()[windingIndex]);
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

    std::vector<std::vector<size_t>> get_patterns(InputsWrapper& inputs) {
        auto isolationSidesRequired = get_isolation_sides(inputs);

        if (!inputs.get_design_requirements().get_isolation_sides()) {
            throw std::runtime_error("Missing isolation sides requirement");
        }

        auto isolationSidesRequirement = inputs.get_design_requirements().get_isolation_sides().value();

        std::vector<std::vector<size_t>> sectionPatterns;
        for(unsigned i = 0; i < tgamma(isolationSidesRequired.size() + 1) / 2; ++i) {
            std::vector<size_t> sectionPattern;
            for (auto isolationSide : isolationSidesRequired) {
                for (size_t windingIndex = 0; windingIndex < inputs.get_mutable_design_requirements().get_turns_ratios().size() + 1; ++windingIndex) {
                    if (isolationSidesRequirement[windingIndex] == isolationSide) {
                        sectionPattern.push_back(windingIndex);
                    }
                }
            }
            sectionPatterns.push_back(sectionPattern);

            std::next_permutation(isolationSidesRequired.begin(), isolationSidesRequired.end());
        }

        return sectionPatterns;
    }


    std::vector<size_t> get_repetitions(InputsWrapper& inputs) {
        if (inputs.get_design_requirements().get_turns_ratios().size() == 0) {
            return {1};  // hardcoded
        }
        return {1, 2};  // hardcoded
        
    }


    std::vector<std::pair<MasWrapper, double>> CoilAdviser::get_advised_coil(MasWrapper mas, size_t maximumNumberResults){
        std::string file_path = __FILE__;
        auto inventory_path = file_path.substr(0, file_path.rfind("/")).append("/../../MAS/data/wires.ndjson");
        std::ifstream ndjsonFile(inventory_path);
        std::string jsonLine;
        std::vector<WireWrapper> wires;
        while (std::getline(ndjsonFile, jsonLine)) {
            json jf = json::parse(jsonLine);
            WireWrapper wire(jf);
            if ((_includeFoil || wire.get_type() != WireType::FOIL) &&
                (_includeRectangular || wire.get_type() != WireType::RECTANGULAR) &&
                (_includeLitz || wire.get_type() != WireType::LITZ) &&
                (_includeRound || wire.get_type() != WireType::ROUND)) {
                wires.push_back(wire);
            }
        }
        return get_advised_coil(&wires, mas, maximumNumberResults);
    }

    std::vector<std::pair<MasWrapper, double>> CoilAdviser::get_advised_coil(std::vector<WireWrapper>* wires, MasWrapper mas, size_t maximumNumberResults){
        auto inputs = mas.get_mutable_inputs();
        auto patterns = get_patterns(inputs);
        auto repetitions = get_repetitions(inputs);
        mas.set_inputs(inputs);

        size_t maximumNumberResultsPerPattern = std::max(2.0, ceil(maximumNumberResults / (patterns.size() * repetitions.size())));

        auto combinationsWithstandVoltageForWires = get_solid_insulation_requirements_for_wires(mas.get_inputs());
        std::vector<std::pair<MasWrapper, double>> masMagneticsWithCoil;
        for (auto repetition : repetitions) {
            for (auto pattern : patterns) {
                for (auto& withstandVoltageForWires : combinationsWithstandVoltageForWires) {
                    auto resultsPerPattern = get_advised_coil_for_pattern(wires, mas, pattern, repetition, withstandVoltageForWires, maximumNumberResultsPerPattern);
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

    std::vector<std::vector<WireSolidInsulationRequirements>> CoilAdviser::get_solid_insulation_requirements_for_wires(InputsWrapper inputs) {
        auto isolationSidesRequired = get_isolation_sides(inputs);
        auto withstandVoltage = InsulationCoordinator().calculate_withstand_voltage(inputs);
        size_t numberWindings = inputs.get_design_requirements().get_turns_ratios().size() + 1;
        auto insulationType = inputs.get_insulation_type();
        bool canFullyInsulatedWireBeUsed = InsulationCoordinator::can_fully_insulated_wire_be_used(inputs);
        std::vector<std::vector<WireSolidInsulationRequirements>> combinationsWithstandVoltageForWires;
        auto isolationSidePerWinding = inputs.get_design_requirements().get_isolation_sides().value();

        // We always add the option of the wires not complying with anything
        {
            std::vector<WireSolidInsulationRequirements> solidInsulationForWires;
            for(unsigned windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
                solidInsulationForWires.push_back(get_requirements_for_functional());
            }
            combinationsWithstandVoltageForWires.push_back(solidInsulationForWires);
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
                    combinationsWithstandVoltageForWires.push_back(solidInsulationForWires);
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
                    for (auto isolationSidePerWinding : isolationSidePerWinding) {
                        solidInsulationForWires.push_back(get_requirements_for_basic(withstandVoltage, canFullyInsulatedWireBeUsed));
                    }
                    combinationsWithstandVoltageForWires.push_back(solidInsulationForWires);
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
                        combinationsWithstandVoltageForWires.push_back(solidInsulationForWires);
                    }
                }
            }
        }

        return combinationsWithstandVoltageForWires;
    }

    std::vector<std::pair<MasWrapper, double>> CoilAdviser::get_advised_coil_for_pattern(std::vector<WireWrapper>* wires, MasWrapper mas, std::vector<size_t> pattern, size_t repetitions, std::vector<WireSolidInsulationRequirements> withstandVoltageForWires, size_t maximumNumberResults){
        auto defaults = Defaults();
        auto sectionProportions = calculate_winding_window_proportion_per_winding(mas.get_inputs());
        size_t numberWindings = mas.get_mutable_magnetic().get_coil().get_functional_description().size();

        // mas.get_mutable_magnetic().get_mutable_coil().set_inputs(mas.get_inputs());
        mas.get_mutable_magnetic().get_mutable_coil().wind_by_sections(sectionProportions, pattern, repetitions);
        mas.get_mutable_magnetic().get_mutable_coil().delimit_and_compact();
        OpenMagnetics::WireAdviser wireAdviser;
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
            wireAdviser.set_wire_solid_insulation_requirements(withstandVoltageForWires[windingIndex]);

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

            for (auto& wireConfiguration : wireConfigurations) {
                wireAdviser.set_maximum_effective_current_density(wireConfiguration["maximumEffectiveCurrentDensity"]);
                wireAdviser.set_maximum_number_parallels(wireConfiguration["maximumNumberParallels"]);

                auto wiresWithScoring = wireAdviser.get_advised_wire(wires,
                                                                     mas.get_mutable_magnetic().get_coil().get_functional_description()[windingIndex],
                                                                     mas.get_mutable_magnetic().get_coil().get_sections_description().value()[windingIndex],
                                                                     maximumCurrent,
                                                                     maximumTemperature,
                                                                     mas.get_mutable_magnetic().get_mutable_coil().get_interleaving_level(),
                                                                     1000);

                if (wiresWithScoring.size() != 0) {
                    timeout += wiresWithScoring.size();
                    wireCoilPerWinding.push_back(wiresWithScoring);
                    break;
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

        auto currentWireIndexPerWinding = std::vector<size_t>(numberWindings, 0);
        std::vector<std::pair<MasWrapper, double>> masMagneticsWithCoil;

        while (true) {

            std::vector<CoilFunctionalDescription> coilFunctionalDescription;

            for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
                coilFunctionalDescription.push_back(wireCoilPerWinding[windingIndex][currentWireIndexPerWinding[windingIndex]].first);
            }
            mas.get_mutable_magnetic().get_mutable_coil().set_functional_description(coilFunctionalDescription);

            bool wound = mas.get_mutable_magnetic().get_mutable_coil().wind(sectionProportions, pattern, repetitions);
            if (wound) {
                mas.get_mutable_magnetic().get_mutable_coil().delimit_and_compact();
                mas.get_mutable_magnetic().set_coil(mas.get_mutable_magnetic().get_mutable_coil());
                masMagneticsWithCoil.push_back({mas, 1.0});
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
