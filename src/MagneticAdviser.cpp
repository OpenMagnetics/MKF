#include "InputsWrapper.h"
#include "MagneticAdviser.h"
#include "MagneticSimulator.h"
#include "Painter.h"
#include "Defaults.h"
#include "Settings.h"
#include <magic_enum_utility.hpp>


namespace OpenMagnetics {

    std::vector<std::pair<MasWrapper, double>> MagneticAdviser::get_advised_magnetic(InputsWrapper inputs, size_t maximumNumberResults) {
        std::map<MagneticAdviserFilters, double> weights;
        magic_enum::enum_for_each<MagneticAdviserFilters>([&] (auto val) {
            MagneticAdviserFilters filter = val;
            weights[filter] = 1.0;
        });
        return get_advised_magnetic(inputs, weights, maximumNumberResults);
    }
    std::vector<std::pair<MasWrapper, double>> MagneticAdviser::get_advised_magnetic(InputsWrapper inputs, std::map<MagneticAdviserFilters, double> weights, size_t maximumNumberResults) {
        auto settings = OpenMagnetics::Settings::GetInstance();
        std::vector<MasWrapper> masData;

        if (coreDatabase.empty()) {
            load_cores(false, settings->get_use_only_cores_in_stock());
        }
        if (wireDatabase.empty()) {
            load_wires();
        }

        std::map<CoreAdviser::CoreAdviserFilters, double> coreWeights;
        coreWeights[CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 1;
        coreWeights[CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 1;
        coreWeights[CoreAdviser::CoreAdviserFilters::COST] = weights[MagneticAdviser::MagneticAdviserFilters::COST];
        coreWeights[CoreAdviser::CoreAdviserFilters::EFFICIENCY] = weights[MagneticAdviser::MagneticAdviserFilters::EFFICIENCY];
        coreWeights[CoreAdviser::CoreAdviserFilters::DIMENSIONS] = weights[MagneticAdviser::MagneticAdviserFilters::DIMENSIONS];

        CoreAdviser coreAdviser(false);
        coreAdviser.set_unique_core_shapes(true);
        CoilAdviser coilAdviser;
        MagneticSimulator magneticSimulator;
        size_t numberWindings = inputs.get_design_requirements().get_turns_ratios().size() + 1;

        double clearanceAndCreepageDistance = InsulationCoordinator().calculate_creepage_distance(inputs, true);
        coreAdviser.set_average_margin_in_winding_window(clearanceAndCreepageDistance);

        // std::cout << "Getting core" << std::endl;
        size_t expectedWoundCores = std::min(maximumNumberResults, std::max(size_t(2), size_t(floor(double(maximumNumberResults) / numberWindings))));
        auto masMagneticsWithCore = coreAdviser.get_advised_core(inputs, coreWeights, expectedWoundCores * 10);
        size_t coresWound = 0;
        for (auto& masMagneticWithCore : masMagneticsWithCore) {
            // std::cout << "core:                                                                 " << masMagneticWithCore.first.get_magnetic().get_core().get_name().value() << std::endl;
            // std::cout << "Getting coil" << std::endl;
            auto masMagneticsWithCoreAndCoil = coilAdviser.get_advised_coil(masMagneticWithCore.first, std::max(2.0, ceil(double(maximumNumberResults) / masMagneticsWithCore.size())));
            if (masMagneticsWithCoreAndCoil.size() > 0) {
                coresWound++;
            }
            size_t processedCoils = 0;
            for (auto mas : masMagneticsWithCoreAndCoil) {

                mas = magneticSimulator.simulate(mas);
                processedCoils++;

                masData.push_back(mas);
                if (processedCoils >= size_t(ceil(maximumNumberResults * 0.5))) {
                    break;
                }
            }
            if (coresWound >= expectedWoundCores) {
                break;
            }
        }

        auto masMagneticsWithScoring = score_magnetics(masData, weights);


        sort(masMagneticsWithScoring.begin(), masMagneticsWithScoring.end(), [](std::pair<MasWrapper, double>& b1, std::pair<MasWrapper, double>& b2) {
            return b1.second > b2.second;
        });

        if (masMagneticsWithScoring.size() > maximumNumberResults) {
            masMagneticsWithScoring = std::vector<std::pair<MasWrapper, double>>(masMagneticsWithScoring.begin(), masMagneticsWithScoring.end() - (masMagneticsWithScoring.size() - maximumNumberResults));
        }

        return masMagneticsWithScoring;

    }

    void MagneticAdviser::normalize_scoring(std::vector<std::pair<MasWrapper, double>>* masMagneticsWithScoring, std::vector<double>* scoring, double weight, std::map<std::string, bool> filterConfiguration) {
        double maximumScoring = *std::max_element(scoring->begin(), scoring->end());
        double minimumScoring = *std::min_element(scoring->begin(), scoring->end());

        for (size_t i = 0; i < (*masMagneticsWithScoring).size(); ++i) {
            auto mas = (*masMagneticsWithScoring)[i].first;
            if (maximumScoring != minimumScoring) {

                if (filterConfiguration["log"]){
                    if (filterConfiguration["invert"]) {
                        (*masMagneticsWithScoring)[i].second += weight * (1 - (std::log10((*scoring)[i]) - std::log10(minimumScoring)) / (std::log10(maximumScoring) - std::log10(minimumScoring)));
                    }
                    else {
                        (*masMagneticsWithScoring)[i].second += weight * (std::log10((*scoring)[i]) - std::log10(minimumScoring)) / (std::log10(maximumScoring) - std::log10(minimumScoring));
                    }
                }
                else {
                    if (filterConfiguration["invert"]) {
                        (*masMagneticsWithScoring)[i].second += weight * (1 - ((*scoring)[i] - minimumScoring) / (maximumScoring - minimumScoring));
                    }
                    else {
                        (*masMagneticsWithScoring)[i].second += weight * ((*scoring)[i] - minimumScoring) / (maximumScoring - minimumScoring);
                    }
                }
            }
            else {
                (*masMagneticsWithScoring)[i].second += 1;
            }
        }
    }

    std::vector<std::pair<MasWrapper, double>> MagneticAdviser::score_magnetics(std::vector<MasWrapper> masMagnetics, std::map<MagneticAdviserFilters, double> weights) {
        _weights = weights;
        std::vector<std::pair<MasWrapper, double>> masMagneticsWithScoring;
        for (auto mas : masMagnetics) {
            masMagneticsWithScoring.push_back({mas, 0.0});
        }
        {
            std::vector<double> scorings;
            for (auto mas : masMagnetics) {
                double scoring = 0;
                for (auto output : mas.get_outputs()) {
                    scoring += output.get_core_losses().value().get_core_losses() + output.get_winding_losses().value().get_winding_losses();
                }
                scorings.push_back(scoring);
                add_scoring(mas.get_magnetic().get_manufacturer_info().value().get_reference().value(), MagneticAdviser::MagneticAdviserFilters::EFFICIENCY, scoring);

            }
            if (masMagneticsWithScoring.size() > 0) {
                normalize_scoring(&masMagneticsWithScoring, &scorings, weights[MagneticAdviserFilters::EFFICIENCY], _filterConfiguration[MagneticAdviserFilters::EFFICIENCY]);
            }
        }

        {
            std::vector<double> scorings;
            for (auto mas : masMagnetics) {
                std::vector<WireWrapper> wires = mas.get_mutable_magnetic().get_mutable_coil().get_wires();
                auto wireRelativeCosts = 0;
                for (auto wire : wires) {
                    wireRelativeCosts += wire.get_relative_cost();
                }
                auto numberLayers = mas.get_mutable_magnetic().get_mutable_coil().get_layers_description()->size();
                double scoring = numberLayers + wireRelativeCosts;


                scorings.push_back(scoring);
                add_scoring(mas.get_magnetic().get_manufacturer_info().value().get_reference().value(), MagneticAdviser::MagneticAdviserFilters::COST, scoring);
            }
            if (masMagneticsWithScoring.size() > 0) {
                normalize_scoring(&masMagneticsWithScoring, &scorings, weights[MagneticAdviserFilters::COST], _filterConfiguration[MagneticAdviserFilters::COST]);
            }
        }
        {
            std::vector<double> scorings;
            for (auto mas : masMagnetics) {
                std::vector<WireWrapper> wires = mas.get_mutable_magnetic().get_mutable_coil().get_wires();
                auto layers = mas.get_mutable_magnetic().get_mutable_coil().get_layers_description().value();
                auto coilDepth = 0.0;
                for (auto layer : layers) {
                    coilDepth += layer.get_dimensions()[0];
                }
                double scoring = coilDepth;

                scorings.push_back(scoring);
                add_scoring(mas.get_magnetic().get_manufacturer_info().value().get_reference().value(), MagneticAdviser::MagneticAdviserFilters::DIMENSIONS, scoring);
            }
            if (masMagneticsWithScoring.size() > 0) {
                normalize_scoring(&masMagneticsWithScoring, &scorings, weights[MagneticAdviserFilters::DIMENSIONS], _filterConfiguration[MagneticAdviserFilters::DIMENSIONS]);
            }
        }
        return masMagneticsWithScoring;
    }


void preview_magnetic(MasWrapper mas) {
    std::string text = "";
    text += "Core shape: " + mas.get_mutable_magnetic().get_mutable_core().get_shape_name() + "\n";
    text += "Core material: " + mas.get_mutable_magnetic().get_mutable_core().get_material_name() + "\n";
    text += "Core gap: " + std::to_string(mas.get_mutable_magnetic().get_mutable_core().get_functional_description().get_gapping()[0].get_length()) + "\n";
    for (size_t windingIndex = 0; windingIndex < mas.get_mutable_magnetic().get_mutable_coil().get_functional_description().size(); ++windingIndex) {
        auto winding = mas.get_mutable_magnetic().get_mutable_coil().get_functional_description()[windingIndex];
        auto wire = OpenMagnetics::CoilWrapper::resolve_wire(winding);
        text += "Winding: " + winding.get_name() + "\n";
        text += "\tNumber Turns: " + std::to_string(winding.get_number_turns()) + "\n";
        text += "\tNumber Parallels: " + std::to_string(winding.get_number_parallels()) + "\n";
        text += "\tWire: " + std::string(magic_enum::enum_name(wire.get_type())) + " " + std::string(magic_enum::enum_name(wire.get_standard().value())) + " " + wire.get_name().value() + "\n";
    }

    for (size_t operatingPointIndex = 0; operatingPointIndex < mas.get_outputs().size(); ++operatingPointIndex) {
        auto output = mas.get_outputs()[operatingPointIndex];
        text += "Operating Point: " + std::to_string(operatingPointIndex + 1) + "\n";
        text += "\tCore losses: " + std::to_string(output.get_core_losses().value().get_core_losses()) + "\n";
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


            text += "\t\t\tOhmic losses: " + std::to_string(windingLosses.get_ohmic_losses().value().get_losses()) + "\n";
            text += "\t\t\tSkin effect losses: " + std::to_string(skinEffectLosses) + "\n";
            text += "\t\t\tProximity effect losses: " + std::to_string(proximityEffectLosses) + "\n";
        }
    }
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

        OpenMagnetics::InputsWrapper inputs(masJson["inputs"]);

        auto outputFolder = inputFilepath.parent_path();


        OpenMagnetics::MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, numberMagnetics);
        for (size_t i = 0; i < masMagnetics.size(); ++i){
            preview_magnetic(masMagnetics[i].first);

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
            auto settings = OpenMagnetics::Settings::GetInstance();
            settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::CONTOUR);

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