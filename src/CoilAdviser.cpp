#include "InputsWrapper.h"
#include "WireAdviser.h"
#include "CoilAdviser.h"


namespace OpenMagnetics {

    std::vector<double> calculate_winding_window_proportion_per_winding(Inputs inputs) {
        auto averagePowerPerWinding = std::vector<double>(inputs.get_operating_points()[0].get_excitations_per_winding().size(), 0);

        for (size_t operationPointIndex = 0; operationPointIndex < inputs.get_operating_points().size(); ++operationPointIndex) {
            for (size_t windingIndex = 0; windingIndex < inputs.get_operating_points()[operationPointIndex].get_excitations_per_winding().size(); ++windingIndex) {
                auto power = InputsWrapper::calculate_instantaneous_power(inputs.get_operating_points()[operationPointIndex].get_excitations_per_winding()[windingIndex]);
                averagePowerPerWinding[windingIndex] += power;
            }
        }

        double sumValue = std::reduce(averagePowerPerWinding.begin(), averagePowerPerWinding.end());
        for (size_t windingIndex = 0; windingIndex < averagePowerPerWinding.size(); ++windingIndex) {
            averagePowerPerWinding[windingIndex] = averagePowerPerWinding[windingIndex] / sumValue;
        }

        return averagePowerPerWinding;
    }

    std::vector<std::pair<MasWrapper, double>> CoilAdviser::get_advised_coil(MasWrapper mas, size_t maximumNumberResults){
        std::cout << "mas.get_inputs().get_operating_points()[0].get_excitations_per_winding().size(): " << mas.get_inputs().get_operating_points()[0].get_excitations_per_winding().size() << std::endl;
        auto sectionProportions = calculate_winding_window_proportion_per_winding(mas.get_inputs());
        std::cout << "sectionProportions.size(): "<< sectionProportions.size() << std::endl;
        std::cout << "sectionProportions[0]: "<< sectionProportions[0] << std::endl;
        std::cout << "sectionProportions[1]: "<< sectionProportions[1] << std::endl;
        mas.get_mutable_magnetic().get_mutable_coil().set_interleaving_level(2);
        mas.get_mutable_magnetic().get_mutable_coil().wind_by_sections(sectionProportions);
        mas.get_mutable_magnetic().get_mutable_coil().delimit_and_compact();

        size_t windingIndex = 0;
        OpenMagnetics::WireAdviser wireAdviser(false, false);

        if (!mas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()) {
            throw std::runtime_error("Missing current in excitaton");
        }

        if (!mas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current().value().get_harmonics() &&
            !mas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current().value().get_processed() &&
            !mas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current().value().get_waveform()) {
            throw std::runtime_error("Missing current harmonics, waveform and processed in excitaton");
        }

        auto coilFunctionalDescriptions = wireAdviser.get_advised_wire(mas.get_mutable_magnetic().get_coil().get_functional_description()[windingIndex],
                                                                       mas.get_mutable_magnetic().get_coil().get_sections_description().value()[windingIndex],
                                                                       mas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current().value(),  // harcoded, should take all op points
                                                                       mas.get_inputs().get_operating_points()[0].get_conditions().get_ambient_temperature(),   // harcoded, should be core or magnetic temperature??
                                                                       mas.get_mutable_magnetic().get_mutable_coil().get_interleaving_level(),
                                                                       1);

        std::cout << "coilFunctionalDescriptions.size(): "<< coilFunctionalDescriptions.size() << std::endl;
        std::cout << "coilFunctionalDescriptions[0].first.get_number_turns(): "<< coilFunctionalDescriptions[0].first.get_number_turns() << std::endl;
        std::cout << "coilFunctionalDescriptions[0].first.get_number_parallels(): "<< coilFunctionalDescriptions[0].first.get_number_parallels() << std::endl;
        std::cout << "CoilWrapper::resolve_wire(coilFunctionalDescriptions[0].first).get_name().value(): "<< CoilWrapper::resolve_wire(coilFunctionalDescriptions[0].first).get_name().value() << std::endl;
        std::cout << "coilFunctionalDescriptions[0].second: "<< coilFunctionalDescriptions[0].second << std::endl;
        std::cout << "CoilWrapper::resolve_wire(coilFunctionalDescriptions[1].first).get_name().value(): "<< CoilWrapper::resolve_wire(coilFunctionalDescriptions[1].first).get_name().value() << std::endl;
        std::cout << "coilFunctionalDescriptions[1].second: "<< coilFunctionalDescriptions[1].second << std::endl;
        std::cout << "CoilWrapper::resolve_wire(coilFunctionalDescriptions[2].first).get_name().value(): "<< CoilWrapper::resolve_wire(coilFunctionalDescriptions[2].first).get_name().value() << std::endl;
        std::cout << "coilFunctionalDescriptions[2].second: "<< coilFunctionalDescriptions[2].second << std::endl;
        std::cout << "CoilWrapper::resolve_wire(coilFunctionalDescriptions[123].first).get_name().value(): "<< CoilWrapper::resolve_wire(coilFunctionalDescriptions[123].first).get_name().value() << std::endl;
        std::cout << "coilFunctionalDescriptions[123].first.get_number_parallels(): "<< coilFunctionalDescriptions[123].first.get_number_parallels() << std::endl;
        // for (auto& ea : coilFunctionalDescriptions) {
        //     std::cout << ea.first.get_number_parallels() << std::endl;

        // }
        return std::vector<std::pair<MasWrapper, double>>{{mas, 1.0}};

    }
} // namespace OpenMagnetics
