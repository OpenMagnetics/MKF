#include "processors/Inputs.h"
#include "advisers/CoilAdviser.h"
#include "Models.h"
#include "constructive_models/Insulation.h"

namespace OpenMagnetics {

    std::vector<double> calculate_winding_window_proportion_per_winding(Inputs& inputs) {
        auto averagePowerPerWinding = std::vector<double>(inputs.get_operating_points()[0].get_excitations_per_winding().size(), 0);

        for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
            for (size_t windingIndex = 0; windingIndex < inputs.get_operating_points()[operatingPointIndex].get_excitations_per_winding().size(); ++windingIndex) {

                if (!inputs.get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_voltage()) {
                    auto excitation = inputs.get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex];
                    inputs.get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex].set_voltage(Inputs::calculate_induced_voltage(excitation, resolve_dimensional_values(inputs.get_design_requirements().get_magnetizing_inductance())));
                }

                auto power = Inputs::calculate_instantaneous_power(inputs.get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex]);
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

    std::vector<Mas> CoilAdviser::get_advised_coil(Mas mas, size_t maximumNumberResults){
        logEntry("Starting Coil Adviser without wires", "CoilAdviser");
        if (wireDatabase.empty()) {
            load_wires();
        }
        std::string jsonLine;
        std::vector<Wire> wires;
        for (const auto& [key, wire] : wireDatabase) {
            if ((settings->get_wire_adviser_include_planar() || wire.get_type() != WireType::PLANAR) &&
                (settings->get_wire_adviser_include_foil() || wire.get_type() != WireType::FOIL) &&
                (settings->get_wire_adviser_include_rectangular() || wire.get_type() != WireType::RECTANGULAR) &&
                (settings->get_wire_adviser_include_litz() || wire.get_type() != WireType::LITZ) &&
                (settings->get_wire_adviser_include_round() || wire.get_type() != WireType::ROUND)) {

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
        auto patterns = Coil::get_patterns(inputs, coreType);
        auto repetitions = Coil::get_repetitions(inputs, coreType);
        mas.set_inputs(inputs);

        size_t maximumNumberResultsPerPattern = std::max(2.0, ceil(maximumNumberResults / (patterns.size() * repetitions.size())));
        logEntry("Trying " + std::to_string(repetitions.size()) + " repetitions and " + std::to_string(patterns.size()) + " patterns", "CoilAdviser");

        std::vector<Mas> masMagneticsWithCoil;
        for (auto repetition : repetitions) {
            for (auto pattern : patterns) {
                auto aux = mas.get_mutable_magnetic().get_mutable_coil().check_pattern_and_repetitions_integrity(pattern, repetition);
                pattern = aux.first;
                repetition = aux.second;
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

                    std::move(resultsPerPattern.begin(), resultsPerPattern.end(), std::back_inserter(masMagneticsWithCoil));
                }
            }
        }

        return masMagneticsWithCoil;
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
            auto sections = coil.get_sections_description().value();
            return sections;
        }
        else {
            return {};
        }
    }

    std::vector<Mas> CoilAdviser::get_advised_coil_for_pattern(std::vector<Wire>* wires, Mas mas, std::vector<size_t> pattern, size_t repetitions, std::vector<WireSolidInsulationRequirements> solidInsulationRequirementsForWires, size_t maximumNumberResults, std::string reference){
        bool filterMode = bool(mas.get_mutable_inputs().get_design_requirements().get_minimum_impedance());
        size_t maximumNumberWires = settings->get_coil_adviser_maximum_number_wires();
        auto sectionProportions = calculate_winding_window_proportion_per_winding(mas.get_mutable_inputs());
        auto core = mas.get_magnetic().get_core();
        auto coil = mas.get_magnetic().get_coil();

        if (mas.get_mutable_inputs().get_wiring_technology() == WiringTechnology::PRINTED) {
            coil.set_winding_orientation(WindingOrientation::CONTIGUOUS);
            coil.set_section_alignment(CoilAlignment::CENTERED);
            coil.set_layers_orientation(WindingOrientation::CONTIGUOUS);
            coil.set_turns_alignment(CoilAlignment::SPREAD);
        }
        else {
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
        }
        mas.get_mutable_magnetic().set_coil(coil);

        size_t numberWindings = coil.get_functional_description().size();

        auto needsMargin = InsulationCoordinator::needs_margin(solidInsulationRequirementsForWires, pattern, repetitions);
        coil.set_inputs(mas.get_inputs());
        coil.clear();
        coil.wind_by_sections(sectionProportions, pattern, repetitions);

        auto sections = get_advised_sections(mas, pattern, repetitions);
        if (sections.size() == 0) {
            return {};
        }
        coil.set_sections_description(sections);

        for (auto section : sections) {
            if (section.get_dimensions()[0] < 0) {
                throw std::runtime_error("section.get_dimensions()[0] cannot be negative");
            }
            if (section.get_dimensions()[1] < 0) {
                throw std::runtime_error("section.get_dimensions()[1] cannot be negative: " + std::to_string(section.get_dimensions()[1]));
            }
        }

        _wireAdviser.set_common_wire_standard(_commonWireStandard);  // TODO, rethink this

        if (needsMargin && mas.get_mutable_inputs().get_wiring_technology() == WiringTechnology::WOUND) {
            // If we want to use margin, we set the maximum so the wires chosen will need margin (and be faster)
            for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
                solidInsulationRequirementsForWires[windingIndex].set_maximum_number_layers(1);
                solidInsulationRequirementsForWires[windingIndex].set_maximum_grade(3);
            }
        }

        std::vector<std::vector<std::pair<OpenMagnetics::CoilFunctionalDescription, double>>> wireCoilPerWinding;


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
                if (!mas.get_inputs().get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current()) {
                    throw std::runtime_error("Current is missing");
                }
                if (!mas.get_inputs().get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current()->get_processed()) {
                    throw std::runtime_error("Current is not processed");
                }
                if (!mas.get_inputs().get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current()->get_processed()->get_effective_frequency()) {
                    auto current = mas.get_inputs().get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex].get_current().value();
                    if (current.get_waveform()) {
                        auto processed = Inputs::calculate_processed_data(current.get_waveform().value());
                        current.set_processed(processed);
                        mas.get_mutable_inputs().get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[windingIndex].set_current(current);
                    }
                    else {
                        throw std::runtime_error("Current is not processed");
                    }
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
            if (mas.get_mutable_inputs().get_wiring_technology() == WiringTechnology::PRINTED) {
                bool found = false;
                auto sectionIndex = coil.convert_conduction_section_index_to_global(windingIndex);

                auto wiresWithScoring = _wireAdviser.get_advised_planar_wire(coil.get_functional_description()[windingIndex],
                                                                             sections[sectionIndex],
                                                                             maximumCurrent,
                                                                             maximumTemperature,
                                                                             coil.get_interleaving_level(),
                                                                             maximumNumberWires);

                if (wiresWithScoring.size() != 0) {
                    timeout += wiresWithScoring.size();

                    wireCoilPerWinding.push_back(wiresWithScoring);
                    found = true;
                }
                if (!found) {
                    wireCoilPerWinding.push_back(std::vector<std::pair<CoilFunctionalDescription, double>>{});
                }
            }
            else {
                std::vector<std::map<std::string, double>> wireConfigurations = {
                    {{"maximumEffectiveCurrentDensity", defaults.maximumEffectiveCurrentDensity}, {"maximumNumberParallels", defaults.maximumNumberParallels}},
                    {{"maximumEffectiveCurrentDensity", defaults.maximumEffectiveCurrentDensity}, {"maximumNumberParallels", defaults.maximumNumberParallels * 2}},
                    {{"maximumEffectiveCurrentDensity", defaults.maximumEffectiveCurrentDensity * 2}, {"maximumNumberParallels", defaults.maximumNumberParallels}},
                    {{"maximumEffectiveCurrentDensity", defaults.maximumEffectiveCurrentDensity * 2}, {"maximumNumberParallels", defaults.maximumNumberParallels * 2}}
                };
                bool found = false;
                logEntry("Trying " + std::to_string(wireConfigurations.size()) + " wire configurations", "CoilAdviser", 2);

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
        }

        logEntry("Trying to wind " + std::to_string(wireCoilPerWinding[0].size()) + " coil possibilities", "CoilAdviser");
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
        std::vector<Mas> masMagneticsWithCoil;

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
        logEntry("Managed to wind " + std::to_string(masMagneticsWithCoil.size()) + " coils", "CoilAdviser");

        return masMagneticsWithCoil;

    }
} // namespace OpenMagnetics
