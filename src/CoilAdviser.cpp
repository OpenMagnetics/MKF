#include "InputsWrapper.h"
#include "WireAdviser.h"
#include "CoilAdviser.h"
#include "Defaults.h"
#include "../tests/TestingUtils.h"
#include "Models.h"
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
        auto defaults = Defaults();
        auto sectionProportions = calculate_winding_window_proportion_per_winding(mas.get_inputs());
        size_t numberWindings = mas.get_mutable_magnetic().get_coil().get_functional_description().size();
        mas.get_mutable_magnetic().get_mutable_coil().set_interleaving_level(_interleavingLevel);
        mas.get_mutable_magnetic().get_mutable_coil().wind_by_sections(sectionProportions);
        {
            std::string filePath = __FILE__;
            auto outputFilePath = filePath.substr(0, filePath.rfind("/")).append("/../output/");
            auto outFile = outputFilePath;
            std::string filename = "Debug_0.svg";
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_core(mas.get_magnetic());
            painter.paint_bobbin(mas.get_magnetic());
            painter.paint_coil_sections(mas.get_magnetic());
            painter.export_svg();
        }
        mas.get_mutable_magnetic().get_mutable_coil().delimit_and_compact();
        {
            std::string filePath = __FILE__;
            auto outputFilePath = filePath.substr(0, filePath.rfind("/")).append("/../output/");
            auto outFile = outputFilePath;
            std::string filename = "Debug_1.svg";
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);
            painter.paint_core(mas.get_magnetic());
            painter.paint_bobbin(mas.get_magnetic());
            painter.paint_coil_sections(mas.get_magnetic());
            painter.export_svg();
        }

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

            // std::cout << "windingIndex: "<< windingIndex << std::endl;
            // std::cout << "Mierdon 1" << std::endl;
            for (auto& wireConfiguration : wireConfigurations) {
                wireAdviser.set_maximum_effective_current_density(wireConfiguration["maximumEffectiveCurrentDensity"]);
                wireAdviser.set_maximum_number_parallels(wireConfiguration["maximumNumberParallels"]);
                // std::cout << "maximumEffectiveCurrentDensity: "<< wireConfiguration["maximumEffectiveCurrentDensity"] << std::endl;
                // std::cout << "maximumNumberParallels: "<< wireConfiguration["maximumNumberParallels"] << std::endl;
                // json mierda;
                // to_json(mierda, mas.get_mutable_magnetic().get_coil().get_functional_description()[windingIndex]);
                // std::cout << "get_functional_description" << std::endl;
                // std::cout << mierda << std::endl;
                // to_json(mierda, mas.get_mutable_magnetic().get_coil().get_sections_description().value()[windingIndex]);
                // std::cout << "get_sections_description" << std::endl;
                // std::cout << mierda << std::endl;

                auto coilsWithScoring = wireAdviser.get_advised_wire(mas.get_mutable_magnetic().get_coil().get_functional_description()[windingIndex],
                                                                               mas.get_mutable_magnetic().get_coil().get_sections_description().value()[windingIndex],
                                                                               maximumCurrent,
                                                                               maximumTemperature,
                                                                               mas.get_mutable_magnetic().get_mutable_coil().get_interleaving_level(),
                                                                               1000);

            // std::cout << "Mierdon 2" << std::endl;
            // std::cout << "coilsWithScoring.size(): " << coilsWithScoring.size() << std::endl;
            // std::cout << "wireAdviser.get_maximum_area_proportion(): " << wireAdviser.get_maximum_area_proportion() << std::endl;
            // std::cout << "sectionProportions[windingIndex]: " << sectionProportions[windingIndex] << std::endl;
                if (coilsWithScoring.size() != 0) {
            // std::cout << "Mierdon 3" << std::endl;
            // std::cout << "maximumNumberParallels: "<< coilsWithScoring[0][0] << std::endl;
                    timeout += coilsWithScoring.size();
                    wireCoilPerWinding.push_back(coilsWithScoring);
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
        bool wound = false;

        while (!wound) {

            std::vector<CoilFunctionalDescription> coilFunctionalDescription;

            for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
                coilFunctionalDescription.push_back(wireCoilPerWinding[windingIndex][currentWireIndexPerWinding[windingIndex]].first);
            }
            mas.get_mutable_magnetic().get_mutable_coil().set_functional_description(coilFunctionalDescription);

            wound = mas.get_mutable_magnetic().get_mutable_coil().wind();
            if (wound) {
                mas.get_mutable_magnetic().get_mutable_coil().delimit_and_compact();
                mas.get_mutable_magnetic().set_coil(mas.get_mutable_magnetic().get_mutable_coil());
                break;
            }
            timeout--;
            if (timeout == 0) {
                // std::cout << "timeout!!!!!!" << std::endl;
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

        std::vector<std::pair<MasWrapper, double>> masMagneticsWithCoil; 

        if (wound) {
            return {{mas, 1.0}};
        }

        return {};

    }
} // namespace OpenMagnetics
