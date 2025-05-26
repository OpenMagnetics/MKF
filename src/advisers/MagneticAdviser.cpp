#include "processors/Inputs.h"
#include "advisers/MagneticAdviser.h"
#include "physical_models/Impedance.h"
#include "processors/MagneticSimulator.h"
#include "support/Painter.h"
#include <magic_enum_utility.hpp>


namespace OpenMagnetics {

std::vector<std::pair<Mas, double>> MagneticAdviser::get_advised_magnetic(Inputs inputs, size_t maximumNumberResults) {
    std::map<MagneticAdviserFilters, double> weights;
    magic_enum::enum_for_each<MagneticAdviserFilters>([&] (auto val) {
        MagneticAdviserFilters filter = val;
        weights[filter] = 1.0;
    });
    return get_advised_magnetic(inputs, weights, maximumNumberResults);
}

std::vector<std::pair<Mas, double>> MagneticAdviser::get_advised_magnetic(Inputs inputs, std::map<MagneticAdviserFilters, double> weights, size_t maximumNumberResults) {
    bool filterMode = bool(inputs.get_design_requirements().get_minimum_impedance());
    std::vector<Mas> masData;

    if (filterMode) {
        settings->set_use_toroidal_cores(true);
        settings->set_use_only_cores_in_stock(false);
        settings->set_use_concentric_cores(false);
    }

    if (coreDatabase.empty()) {
        load_cores();
    }
    if (wireDatabase.empty()) {
        load_wires();
    }

    bool previousCoilIncludeAdditionalCoordinates = settings->get_coil_include_additional_coordinates();
    settings->set_coil_include_additional_coordinates(false);

    std::map<CoreAdviser::CoreAdviserFilters, double> coreWeights;
    coreWeights[CoreAdviser::CoreAdviserFilters::COST] = weights[MagneticAdviser::MagneticAdviserFilters::COST];
    coreWeights[CoreAdviser::CoreAdviserFilters::DIMENSIONS] = weights[MagneticAdviser::MagneticAdviserFilters::DIMENSIONS];

    if (filterMode) {
        coreWeights[CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 0;
        coreWeights[CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 0;
        coreWeights[CoreAdviser::CoreAdviserFilters::MINIMUM_IMPEDANCE] = weights[MagneticAdviser::MagneticAdviserFilters::EFFICIENCY];
        coreWeights[CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 0;
    }
    else {
        coreWeights[CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 1;
        coreWeights[CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 1;
        coreWeights[CoreAdviser::CoreAdviserFilters::EFFICIENCY] = weights[MagneticAdviser::MagneticAdviserFilters::EFFICIENCY];
        coreWeights[CoreAdviser::CoreAdviserFilters::MINIMUM_IMPEDANCE] = 0;
    }

    CoreAdviser coreAdviser;
    coreAdviser.set_unique_core_shapes(true);
    CoilAdviser coilAdviser;
    MagneticSimulator magneticSimulator;
    size_t numberWindings = inputs.get_design_requirements().get_turns_ratios().size() + 1;

    std::cout << "Getting core" << std::endl;
    size_t expectedWoundCores = std::min(maximumNumberResults, std::max(size_t(2), size_t(floor(double(maximumNumberResults) / numberWindings))));
    auto masMagneticsWithCore = coreAdviser.get_advised_core(inputs, coreWeights, expectedWoundCores * 10);
    
    size_t coresWound = 0;
    for (auto& [core, coreScoring] : masMagneticsWithCore) {
        std::cout << "core:                                                                 " << core.get_magnetic().get_core().get_name().value() << std::endl;
        std::cout << "Getting coil" << std::endl;
        std::vector<std::pair<size_t, double>> usedNumberSectionsAndMargin;

        auto masMagneticsWithCoreAndCoil = coilAdviser.get_advised_coil(core, std::max(2.0, ceil(double(maximumNumberResults) / masMagneticsWithCore.size())));
        if (masMagneticsWithCoreAndCoil.size() > 0) {
            std::cout << "Core wound!" << std::endl;

            coresWound++;
        }
        size_t processedCoils = 0;
        for (auto mas : masMagneticsWithCoreAndCoil) {
            size_t numberSections = mas.get_magnetic().get_coil().get_sections_description()->size();
            double margin = mas.get_magnetic().get_coil().get_sections_description().value()[0].get_margin().value()[0];
            std::pair<size_t, double> numberSectionsAndMarginCombination = {numberSections, margin};
            if (std::find(usedNumberSectionsAndMargin.begin(), usedNumberSectionsAndMargin.end(), numberSectionsAndMarginCombination) != usedNumberSectionsAndMargin.end()) {
                continue;
            }

            if (previousCoilIncludeAdditionalCoordinates) {
                settings->set_coil_include_additional_coordinates(previousCoilIncludeAdditionalCoordinates);
                mas.get_mutable_magnetic().get_mutable_coil().delimit_and_compact();
                settings->set_coil_include_additional_coordinates(false);
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


    auto masMagneticsWithScoring = score_magnetics(masData, weights);

    sort(masMagneticsWithScoring.begin(), masMagneticsWithScoring.end(), [](std::pair<Mas, double>& b1, std::pair<Mas, double>& b2) {
        return b1.second > b2.second;
    });

    if (masMagneticsWithScoring.size() > maximumNumberResults) {
        masMagneticsWithScoring = std::vector<std::pair<Mas, double>>(masMagneticsWithScoring.begin(), masMagneticsWithScoring.end() - (masMagneticsWithScoring.size() - maximumNumberResults));
    }

    settings->set_coil_include_additional_coordinates(previousCoilIncludeAdditionalCoordinates);
    return masMagneticsWithScoring;
}

std::vector<std::pair<Mas, double>> MagneticAdviser::get_advised_magnetic(Inputs inputs, std::vector<Magnetic> catalogMagnetics, size_t maximumNumberResults) {
    std::vector<Mas> validMagnetics;
    MagneticSimulator magneticSimulator;

    std::map<std::string, std::map<std::string, double>> scoringPerReferencePerRequirement;
    std::map<std::string, double> scoringPerReference;
    std::vector<Magnetic> catalogMagneticsWithSameTurnsRatio;
    MagneticFilterTurnsRatios _magneticFilterTurnsRatios;
    MagneticFilterMaximumDimensions _magneticFilterMaximumDimensions;
    MagneticFilterSaturation _magneticFilterSaturation;
    MagneticFilterDcCurrentDensity _magneticFilterDcCurrentDensity;
    MagneticFilterEffectiveCurrentDensity _magneticFilterEffectiveCurrentDensity;
    MagneticFilterImpedance _magneticFilterImpedance;
    MagneticFilterMagnetizingInductance _magneticFilterMagnetizingInductance;

    for (size_t magneticIndex = 0; magneticIndex < catalogMagnetics.size(); ++magneticIndex) {
        auto magnetic = catalogMagnetics[magneticIndex];
        std::string reference = magnetic.get_manufacturer_info().value().get_reference().value();
        scoringPerReference[reference] = 0;

        // Turns ratios
        {
            auto [valid, scoring] = _magneticFilterTurnsRatios.evaluate_magnetic(&magnetic, &inputs);
            if (!valid)
                continue;
        }

        catalogMagneticsWithSameTurnsRatio.push_back(magnetic);

        // Maximum dimensions
        {
            auto [valid, scoring] = _magneticFilterMaximumDimensions.evaluate_magnetic(&magnetic, &inputs);
            scoringPerReferencePerRequirement["maximumDimensions"][reference] = scoring;
            if (!valid)
                continue;
        }

        // Saturation
        {
            auto [valid, scoring] = _magneticFilterSaturation.evaluate_magnetic(&magnetic, &inputs);
            scoringPerReferencePerRequirement["saturation"][reference] = scoring;
            if (!valid)
                continue;
        }

        // DC Current Density
        {
            auto [valid, scoring] = _magneticFilterDcCurrentDensity.evaluate_magnetic(&magnetic, &inputs);
            scoringPerReferencePerRequirement["dcCurrentDensity"][reference] = scoring;
            if (!valid)
                continue;
        }

        // Effective Current Density
        {
            auto [valid, scoring] = _magneticFilterEffectiveCurrentDensity.evaluate_magnetic(&magnetic, &inputs);
            scoringPerReferencePerRequirement["effectiveCurrentDensity"][reference] = scoring;
            if (!valid)
                continue;
        }

        // if (inputs.get_design_requirements().get_maximum_weight()) {
            // Nice to have in the future
        // }


        // Impedance
        {
            auto [valid, scoring] = _magneticFilterImpedance.evaluate_magnetic(&magnetic, &inputs);
            scoringPerReferencePerRequirement["impedance"][reference] = scoring;
            if (!valid)
                continue;
        }

        // Magnetizing Inductance
        {
            auto [valid, scoring] = _magneticFilterMagnetizingInductance.evaluate_magnetic(&magnetic, &inputs);
            scoringPerReferencePerRequirement["magnetizingInductance"][reference] = scoring;
            if (!valid)
                continue;
        }


        // if (inputs.get_design_requirements().get_operating_temperature()) {
            // Nice to have in the future
        // }
        // if (inputs.get_design_requirements().get_stray_capacitance()) {
            // Nice to have in the future
        // }
        // if (inputs.get_design_requirements().get_insulation()) {
            // Nice to have in the future
        // }

        Mas mas;
        mas.set_magnetic(magnetic);
        mas.set_inputs(inputs);
        mas = magneticSimulator.simulate(mas, true);
        validMagnetics.push_back(mas);
    }

    std::vector<std::pair<Mas, double>> masMagneticsWithScoring;

    if (validMagnetics.size() > 0) {
        std::map<MagneticAdviserFilters, double> weights;
        magic_enum::enum_for_each<MagneticAdviserFilters>([&] (auto val) {
            MagneticAdviserFilters filter = val;
            weights[filter] = 1.0;
        });

        masMagneticsWithScoring = score_magnetics(validMagnetics, weights);
        sort(masMagneticsWithScoring.begin(), masMagneticsWithScoring.end(), [](std::pair<Mas, double>& b1, std::pair<Mas, double>& b2) {
            return b1.second > b2.second;
        });
        if (masMagneticsWithScoring.size() > maximumNumberResults) {
            masMagneticsWithScoring = std::vector<std::pair<Mas, double>>(masMagneticsWithScoring.begin(), masMagneticsWithScoring.end() - (masMagneticsWithScoring.size() - maximumNumberResults));
        }

        return masMagneticsWithScoring;

    }
    else {
        for (auto [filter, aux] : scoringPerReferencePerRequirement) {
            auto normalizedScoring = OpenMagnetics::normalize_scoring(aux, 1, { {"invert", true}, {"log", false} });
            for (auto [reference, scoring] : normalizedScoring) {
                scoringPerReference[reference] += scoring;
            }
        }

        std::vector<std::pair<Mas, double>> masMagneticsWithScoringNoSimulation;
        for (auto magnetic : catalogMagneticsWithSameTurnsRatio) {
            std::string reference = magnetic.get_manufacturer_info().value().get_reference().value();
            Mas mas;
            mas.set_magnetic(magnetic);
            mas.set_inputs(inputs);
            masMagneticsWithScoringNoSimulation.push_back({mas, -scoringPerReference[reference]});
        }
        sort(masMagneticsWithScoringNoSimulation.begin(), masMagneticsWithScoringNoSimulation.end(), [](std::pair<Mas, double>& b1, std::pair<Mas, double>& b2) {
            return b1.second < b2.second;
        });

        if (masMagneticsWithScoringNoSimulation.size() > maximumNumberResults) {
            masMagneticsWithScoringNoSimulation = std::vector<std::pair<Mas, double>>(masMagneticsWithScoringNoSimulation.begin(), masMagneticsWithScoringNoSimulation.end() - (masMagneticsWithScoringNoSimulation.size() - maximumNumberResults));
        }

        for (auto [mas, scoring] : masMagneticsWithScoringNoSimulation) {
            mas = magneticSimulator.simulate(mas, true);
            masMagneticsWithScoring.push_back({mas, scoring});
        }

        return masMagneticsWithScoring;
    }
}

void MagneticAdviser::normalize_scoring(std::vector<std::pair<Mas, double>>* masMagneticsWithScoring, std::vector<double> scoring, double weight, std::map<std::string, bool> filterConfiguration) {
    auto normalizedScorings = OpenMagnetics::normalize_scoring(scoring, weight, filterConfiguration);

    for (size_t i = 0; i < (*masMagneticsWithScoring).size(); ++i) {
        (*masMagneticsWithScoring)[i].second += normalizedScorings[i];
    }
}

std::vector<std::pair<Mas, double>> MagneticAdviser::score_magnetics(std::vector<Mas> masMagnetics, std::map<MagneticAdviserFilters, double> weights) {
    _weights = weights;
    std::vector<std::pair<Mas, double>> masMagneticsWithScoring;
    for (auto mas : masMagnetics) {
        masMagneticsWithScoring.push_back({mas, 0.0});
    }


    MagneticFilterLosses _magneticFilterLosses;
    MagneticFilterCost _magneticFilterCost;
    MagneticFilterDimensions _magneticFilterDimensions;

    {
        std::vector<double> scorings;
        for (auto mas : masMagnetics) {
            auto [valid, scoring] = _magneticFilterLosses.evaluate_magnetic(&mas.get_mutable_magnetic(), &mas.get_mutable_inputs());

            scorings.push_back(scoring);
            add_scoring(mas.get_magnetic().get_manufacturer_info().value().get_reference().value(), MagneticAdviser::MagneticAdviserFilters::EFFICIENCY, scoring);

        }
        if (masMagneticsWithScoring.size() > 0) {
            normalize_scoring(&masMagneticsWithScoring, scorings, weights[MagneticAdviserFilters::EFFICIENCY], _filterConfiguration[MagneticAdviserFilters::EFFICIENCY]);
        }
    }

    {
        std::vector<double> scorings;
        for (auto mas : masMagnetics) {
            auto [valid, scoring] = _magneticFilterCost.evaluate_magnetic(&mas.get_mutable_magnetic(), &mas.get_mutable_inputs());

            scorings.push_back(scoring);
            add_scoring(mas.get_magnetic().get_manufacturer_info().value().get_reference().value(), MagneticAdviser::MagneticAdviserFilters::COST, scoring);
        }
        if (masMagneticsWithScoring.size() > 0) {
            normalize_scoring(&masMagneticsWithScoring, scorings, weights[MagneticAdviserFilters::COST], _filterConfiguration[MagneticAdviserFilters::COST]);
        }
    }
    {
        std::vector<double> scorings;
        for (auto mas : masMagnetics) {
            auto [valid, scoring] = _magneticFilterDimensions.evaluate_magnetic(&mas.get_mutable_magnetic(), &mas.get_mutable_inputs());

            scorings.push_back(scoring);
            add_scoring(mas.get_magnetic().get_manufacturer_info().value().get_reference().value(), MagneticAdviser::MagneticAdviserFilters::DIMENSIONS, scoring);
        }
        if (masMagneticsWithScoring.size() > 0) {
            normalize_scoring(&masMagneticsWithScoring, scorings, weights[MagneticAdviserFilters::DIMENSIONS], _filterConfiguration[MagneticAdviserFilters::DIMENSIONS]);
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
    text += "Core stacks: " + std::to_string(mas.get_mutable_magnetic().get_mutable_core().get_functional_description().get_number_stacks().value()) + "\n";
    for (size_t windingIndex = 0; windingIndex < mas.get_mutable_magnetic().get_mutable_coil().get_functional_description().size(); ++windingIndex) {
        auto winding = mas.get_mutable_magnetic().get_mutable_coil().get_functional_description()[windingIndex];
        auto wire = OpenMagnetics::Coil::resolve_wire(winding);
        text += "Winding: " + winding.get_name() + "\n";
        text += "\tNumber Turns: " + std::to_string(winding.get_number_turns()) + "\n";
        text += "\tNumber Parallels: " + std::to_string(winding.get_number_parallels()) + "\n";
        text += "\tWire: " + std::string(magic_enum::enum_name(wire.get_type()));
        if (wire.get_standard()) {
            text += " " + std::string(magic_enum::enum_name(wire.get_standard().value()));
        }
        if (wire.get_name()) {
            text += " " + wire.get_name().value();
        }
        text += "\n";
    }

    for (size_t operatingPointIndex = 0; operatingPointIndex < mas.get_outputs().size(); ++operatingPointIndex) {
        auto output = mas.get_outputs()[operatingPointIndex];
        text += "Operating Point: " + std::to_string(operatingPointIndex + 1) + "\n";
        text += "\tMagnetizing Inductance: " + std::to_string(resolve_dimensional_values(output.get_magnetizing_inductance().value().get_magnetizing_inductance())) + "\n";
        text += "\tCore losses: " + std::to_string(output.get_core_losses().value().get_core_losses()) + "\n";
        text += "\tMagnetic flux density: " + std::to_string(output.get_core_losses().value().get_magnetic_flux_density()->get_processed()->get_peak().value()) + "\n";
        text += "\tCore temperature: " + std::to_string(output.get_core_losses().value().get_temperature().value()) + "\n";
        text += "\tWinding losses: " + std::to_string(output.get_winding_losses().value().get_winding_losses()) + "\n";
        for (size_t windingIndex = 0; windingIndex < output.get_winding_losses().value().get_winding_losses_per_winding().value().size(); ++windingIndex) {
            auto winding = mas.get_mutable_magnetic().get_mutable_coil().get_functional_description()[windingIndex];
            auto windingLosses = output.get_winding_losses().value().get_winding_losses_per_winding().value()[windingIndex];
            text += "\t\tLosses for winding: " + winding.get_name() + "\n";
            double skinEffectLosses = 0;
            double proximityEffectLosses = 0;
            for (size_t i = 0; i < windingLosses.get_skin_effect_losses().value().get_losses_per_harmonic().size(); ++i)
            {
                skinEffectLosses += windingLosses.get_skin_effect_losses().value().get_losses_per_harmonic()[i];
            }
            for (size_t i = 0; i < windingLosses.get_proximity_effect_losses().value().get_losses_per_harmonic().size(); ++i)
            {
                proximityEffectLosses += windingLosses.get_proximity_effect_losses().value().get_losses_per_harmonic()[i];
            }

            text += "\t\t\tDC resistance: " + std::to_string(output.get_winding_losses().value().get_dc_resistance_per_winding().value()[windingIndex]) + "\n";
            text += "\t\t\tOhmic losses: " + std::to_string(windingLosses.get_ohmic_losses().value().get_losses()) + "\n";
            text += "\t\t\tSkin effect losses: " + std::to_string(skinEffectLosses) + "\n";
            text += "\t\t\tProximity effect losses: " + std::to_string(proximityEffectLosses) + "\n";

            if (windingIndex > 0) {
                std::cout <<  output.get_leakage_inductance().value().get_leakage_inductance_per_winding().size() << std::endl;
                // output.get_leakage_inductance().value().get_leakage_inductance_per_winding();
                double value = output.get_leakage_inductance().value().get_leakage_inductance_per_winding()[windingIndex - 1].get_nominal().value();
                text += "\t\t\tLeakage inductance referred to primary: " + std::to_string(value) + "\n";
            }
        }
    }
    std::cout << text << std::endl;
}

std::map<std::string, std::map<MagneticAdviser::MagneticAdviserFilters, double>> MagneticAdviser::get_scorings(bool weighted){
    std::map<std::string, std::map<MagneticAdviser::MagneticAdviserFilters, double>> swappedScorings;
    for (auto& [filter, aux] : _scorings) {
        auto filterConfiguration = _filterConfiguration[filter];

        double maximumScoring = (*std::max_element(aux.begin(), aux.end(),
                                     [](const std::pair<std::string, double> &p1,
                                        const std::pair<std::string, double> &p2)
                                     {
                                         return p1.second < p2.second;
                                     })).second; 
        double minimumScoring = (*std::min_element(aux.begin(), aux.end(),
                                     [](const std::pair<std::string, double> &p1,
                                        const std::pair<std::string, double> &p2)
                                     {
                                         return p1.second < p2.second;
                                     })).second; 

        for (auto& [name, scoring] : aux) {
            if (filterConfiguration["log"]){
                if (filterConfiguration["invert"]) {
                    if (weighted) {
                        swappedScorings[name][filter] = _weights[filter] * (1 - (std::log10(scoring) - std::log10(minimumScoring)) / (std::log10(maximumScoring) - std::log10(minimumScoring)));
                    }
                    else {
                        swappedScorings[name][filter] = 1 - (std::log10(scoring) - std::log10(minimumScoring)) / (std::log10(maximumScoring) - std::log10(minimumScoring));
                    }
                }
                else {
                    if (weighted) {
                        swappedScorings[name][filter] = _weights[filter] * (std::log10(scoring) - std::log10(minimumScoring)) / (std::log10(maximumScoring) - std::log10(minimumScoring));
                    }
                    else {
                        swappedScorings[name][filter] = (std::log10(scoring) - std::log10(minimumScoring)) / (std::log10(maximumScoring) - std::log10(minimumScoring));
                    }
                }
            }
            else {
                if (filterConfiguration["invert"]) {
                    if (weighted) {
                        swappedScorings[name][filter] = _weights[filter] * (1 - (scoring - minimumScoring) / (maximumScoring - minimumScoring));
                    }
                    else {
                        swappedScorings[name][filter] = 1 - (scoring - minimumScoring) / (maximumScoring - minimumScoring);
                    }
                }
                else {
                    if (weighted) {
                        swappedScorings[name][filter] = _weights[filter] * (scoring - minimumScoring) / (maximumScoring - minimumScoring);
                    }
                    else {
                        swappedScorings[name][filter] = (scoring - minimumScoring) / (maximumScoring - minimumScoring);
                    }
                }
            }
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

int main(int argc, char* argv[]) {
    if (argc == 1) {
        throw std::invalid_argument("Missing inputs file");
    }
    else {
        int numberMagnetics = 1;
        std::filesystem::path inputFilepath = argv[1];
        auto outputFilePath = std::filesystem::path {__FILE__}.parent_path().append("..").append("output");
        if (argc >= 3) {
            if (is_number(argv[2])) {
                numberMagnetics = std::stoi(argv[2]);
            }
            else {
                outputFilePath = argv[2];
            }
        }
        if (argc >= 4) {
            if (!is_number(argv[2]) && is_number(argv[3])) {
                numberMagnetics = std::stoi(argv[3]);
            }
        }

        std::ifstream f(inputFilepath);
        json masJson;
        std::string str;
        if(static_cast<bool>(f)) {
            std::ostringstream ss;
            ss << f.rdbuf(); // reading data
            masJson = json::parse(ss.str());
        }

        if (masJson["inputs"]["designRequirements"]["insulation"]["cti"] == "Group IIIa") {
            masJson["inputs"]["designRequirements"]["insulation"]["cti"] = "Group IIIA";
        }
        if (masJson["inputs"]["designRequirements"]["insulation"]["cti"] == "Group IIIb") {
            masJson["inputs"]["designRequirements"]["insulation"]["cti"] = "Group IIIB";
        }

        OpenMagnetics::Inputs inputs(masJson["inputs"]);

        auto outputFolder = inputFilepath.parent_path();


        OpenMagnetics::MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, numberMagnetics);
        for (size_t i = 0; i < masMagnetics.size(); ++i){
            OpenMagnetics::MagneticAdviser::preview_magnetic(masMagnetics[i].first);

            std::filesystem::path outputFilename = outputFilePath;
            outputFilename += inputFilepath.filename();
            outputFilename += "_design_" + std::to_string(i) + ".json";
            std::ofstream outputFile(outputFilename);
            
            json output;
            to_json(output, masMagnetics[i].first);
            outputFile << output;
            
            outputFile.close();

            outputFilename.replace_extension("svg");
            OpenMagnetics::Painter painter(outputFilename);
            settings->set_painter_mode(OpenMagnetics::PainterModes::CONTOUR);

            settings->set_painter_number_points_x(20);
            settings->set_painter_number_points_y(20);
            settings->set_painter_include_fringing(false);
            settings->set_painter_mirroring_dimension(0);

            painter.paint_magnetic_field(masMagnetics[i].first.get_mutable_inputs().get_operating_point(0), masMagnetics[i].first.get_mutable_magnetic());
            painter.paint_core(masMagnetics[i].first.get_mutable_magnetic());
            painter.paint_bobbin(masMagnetics[i].first.get_mutable_magnetic());
            painter.paint_coil_turns(masMagnetics[i].first.get_mutable_magnetic());
            painter.export_svg();
        }


    }
}