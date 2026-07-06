#include "constructive_models/Wire.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Magnetic.h"
#include "physical_models/MagnetizingInductance.h"
#include "processors/Inputs.h"
#include "support/Utils.h"
#include "advisers/MagneticCombinator.h"
#include <MAS.hpp>
#include <algorithm>
#include <cfloat>
#include <string>

namespace OpenMagnetics {

namespace {

bool repeated_windings(std::vector<size_t> leftWindings, std::vector<size_t> rightWindings) {
    bool repeated = false;
    for (auto windingIndex : leftWindings) {
        if (std::find(rightWindings.begin(), rightWindings.end(), windingIndex) != rightWindings.end()) {
            repeated = true;
            break;
        }
    }

    return repeated;
}

std::string vector_to_string(std::vector<size_t> vector) {
    std::string string;
    for (auto elem : vector) {
        string += std::to_string(elem);
        string += ",";
    }
    string.pop_back();
    return string;
}

std::vector<std::string> get_compulsory_series_connections(size_t windingIndex, std::vector<MagneticCombinator::WindingCombination> combinations) {
    std::vector<std::string> compulsorySeriesConnections;
    auto windingToCheck = combinations[windingIndex];
    if (!windingToCheck.get_connections()) {
        return {};
    }

    auto connectionsToCheck = windingToCheck.get_connections().value();

    for (auto connection : connectionsToCheck) {
        if (connection.get_direction() && connection.get_pin_name()) {
            for (size_t combinationIndex = 0; combinationIndex < combinations.size(); combinationIndex++) {
                if (combinationIndex == windingIndex) {
                    continue;
                }
                if (!combinations[combinationIndex].get_connections()) {
                    continue;
                }
                auto connectionsToCheckAgainst = combinations[combinationIndex].get_connections().value();
                for (auto connectionToCheckAgainst : connectionsToCheckAgainst) {
                    if (connectionToCheckAgainst.get_direction() && connectionToCheckAgainst.get_pin_name()) {
                        if (connectionToCheckAgainst.get_direction().value() != connection.get_direction().value() && connection.get_pin_name().value() == connectionToCheckAgainst.get_pin_name().value()) {
                            compulsorySeriesConnections.push_back(combinations[combinationIndex].get_name());
                        }
                    }
                }
            }
        }
    }

    return compulsorySeriesConnections;

}

} // namespace

OpenMagnetics::Wire MagneticCombinator::get_most_restrictive_wire(OpenMagnetics::Coil coil, std::vector<size_t> windingIndexes) {
    OpenMagnetics::Wire mostRestrictiveWire;
    double minimumDimension = DBL_MAX;
    for (auto windingIndex : windingIndexes) {
        auto wire = coil.resolve_wire(windingIndex);
        double dimension = 0;
        switch (wire.get_type())
        {
            case MAS::WireType::ROUND:
                dimension = OpenMagnetics::resolve_dimensional_values(wire.get_conducting_diameter().value());
                break;
            case MAS::WireType::LITZ:
                {
                    auto strand = wire.resolve_strand();
                    dimension = OpenMagnetics::Wire::get_outer_diameter_served_litz(OpenMagnetics::resolve_dimensional_values(strand.get_conducting_diameter()), wire.get_number_conductors().value());
                    break;
                }
            case MAS::WireType::RECTANGULAR:
                dimension = std::max(OpenMagnetics::resolve_dimensional_values(wire.get_conducting_width().value()), OpenMagnetics::resolve_dimensional_values(wire.get_conducting_height().value()));
                break;
            case MAS::WireType::FOIL:
                dimension = OpenMagnetics::resolve_dimensional_values(wire.get_conducting_width().value());
                break;
        }
        if (dimension < minimumDimension)
        {
            mostRestrictiveWire = wire;
            minimumDimension = dimension;
        }
    }
    return mostRestrictiveWire;
}

std::vector<size_t> MagneticCombinator::get_winding_indexes(OpenMagnetics::Magnetic magnetic) {
    std::vector<size_t> currentWindings;
    for (auto winding : magnetic.get_coil().get_functional_description()) {
        for (auto elem : OpenMagnetics::split(winding.get_name(), ",")) {
            currentWindings.push_back(std::stoul(elem));
        }
    }
    return currentWindings;
}

std::pair<std::map<std::string, MagneticCombinator::WindingCombination>, std::vector<OpenMagnetics::Magnetic>> MagneticCombinator::calculate_virtual_magnetics(OpenMagnetics::Magnetic magnetic, OpenMagnetics::Inputs inputs) {
    auto coil = magnetic.get_coil();
    auto combinations = calculate_coil_combinations(coil, (inputs.get_design_requirements().get_turns_ratios().size() + 1));
    std::vector<OpenMagnetics::Magnetic> virtualMagnetics;
    magnetic.get_mutable_core().process_data();
    magnetic.get_mutable_core().process_gap();

    std::map<std::string, MagneticCombinator::WindingCombination> usedCombinations;

    std::vector<MagneticCombinator::WindingCombination> primaryCombinations;
    for (auto combination : combinations) {
        if (combination.get_isolation_side() == MAS::IsolationSide::PRIMARY) {
            primaryCombinations.push_back(combination);  // TODO: review logic of getting only primaries here
        }
    }

    for (auto combination : primaryCombinations) {
        OpenMagnetics::Magnetic virtualMagnetic;
        virtualMagnetic.set_core(magnetic.get_core());

        OpenMagnetics::Coil virtualCoil;
        virtualCoil.set_bobbin(coil.get_bobbin());
        OpenMagnetics::Winding winding;

        winding.set_connections(combination.get_connections());
        winding.set_isolation_side(combination.get_isolation_side());
        winding.set_number_turns(combination.get_number_turns());
        winding.set_number_parallels(combination.get_number_parallels());
        winding.set_wire(get_most_restrictive_wire(coil, combination.get_winding_indexes()));
        winding.set_name(vector_to_string(combination.get_winding_indexes()));
        virtualCoil.get_mutable_functional_description().push_back(winding);
        virtualMagnetic.set_coil(virtualCoil);

        auto magnetizingInductance = OpenMagnetics::MagnetizingInductance();
        MagnetizingInductanceOutput auxOutput;
        if (inputs.get_operating_points().size() > 0) {
            auto operatingPoint = inputs.get_operating_points()[0];
            auxOutput = magnetizingInductance.calculate_inductance_from_number_turns_and_gapping(virtualMagnetic.get_core(), virtualMagnetic.get_coil(), &operatingPoint);
        }
        else {
            auxOutput = magnetizingInductance.calculate_inductance_from_number_turns_and_gapping(virtualMagnetic.get_core(), virtualMagnetic.get_coil());
        }
        double magnetizingInductanceValue = OpenMagnetics::resolve_dimensional_values(auxOutput.get_magnetizing_inductance());
        bool isValid = OpenMagnetics::check_requirement(inputs.get_design_requirements().get_magnetizing_inductance(), magnetizingInductanceValue);
        if (isValid) {
            virtualMagnetics.push_back(virtualMagnetic);
            usedCombinations[winding.get_name()]  = combination;
        }
    }

    for (auto turnsRatio : inputs.get_design_requirements().get_turns_ratios()) {
        std::vector<OpenMagnetics::Magnetic> temporalVirtualMagnetics;
        for (auto virtualMagnetic : virtualMagnetics) {
            std::vector<size_t> currentWindings = get_winding_indexes(virtualMagnetic);

            for (auto combination : combinations) {
                if (!repeated_windings(currentWindings, combination.get_winding_indexes()))
                {
                    double combinationTurnsRatio = double(virtualMagnetic.get_coil().get_functional_description()[0].get_number_turns()) / combination.get_number_turns();
                    bool isValid = OpenMagnetics::check_requirement(turnsRatio, combinationTurnsRatio);
                    if (isValid) {
                        auto newVirtualMagnetic = virtualMagnetic;
                        OpenMagnetics::Winding winding;
                        winding.set_connections(combination.get_connections());
                        winding.set_isolation_side(combination.get_isolation_side());
                        winding.set_number_turns(combination.get_number_turns());
                        winding.set_number_parallels(combination.get_number_parallels());
                        winding.set_wire(get_most_restrictive_wire(coil, combination.get_winding_indexes()));
                        winding.set_name(vector_to_string(combination.get_winding_indexes()));
                        newVirtualMagnetic.get_mutable_coil().get_mutable_functional_description().push_back(winding);
                        std::vector<size_t> newCurrentWindings = get_winding_indexes(newVirtualMagnetic);
                        bool found = false;
                        for (auto tempMagnetic : temporalVirtualMagnetics) {
                            std::vector<size_t> tempWindingIndexes = get_winding_indexes(tempMagnetic);
                            if (tempWindingIndexes.size() == newCurrentWindings.size()) {
                                bool allIndexesEqual = true;
                                for (size_t index = 0; index < tempWindingIndexes.size(); index++) {
                                    if (tempWindingIndexes[index] != newCurrentWindings[index]) {
                                        allIndexesEqual = false;
                                        break;
                                    }
                                }

                                if (allIndexesEqual) {
                                    found = true;
                                    break;
                                }
                            }

                        }
                        if (!found) {
                            temporalVirtualMagnetics.push_back(newVirtualMagnetic);
                            usedCombinations[winding.get_name()]  = combination;
                        }
                    }
                }
            }
        }

        virtualMagnetics = temporalVirtualMagnetics;
    }

    for (size_t index = 0; index < virtualMagnetics.size(); index++)
    {
        MagneticManufacturerInfo manufacturerInfo;
        if (magnetic.get_manufacturer_info()) {
            manufacturerInfo = magnetic.get_manufacturer_info().value();
        }
        manufacturerInfo.set_name("Wuerth Elektronik");
        manufacturerInfo.set_reference(magnetic.get_manufacturer_info()->get_reference().value() + " virtual magnetic " + std::to_string(index));
        virtualMagnetics[index].set_manufacturer_info(manufacturerInfo);
    }

    return {usedCombinations, virtualMagnetics};
}

OpenMagnetics::Mas MagneticCombinator::devirtualize_mas(OpenMagnetics::Magnetic originalMagnetic, OpenMagnetics::Mas virtualMas, std::map<std::string, MagneticCombinator::WindingCombination> devirtualizationData) {
    auto virtualMagnetic = virtualMas.get_magnetic();
    auto inputs = virtualMas.get_inputs();
    auto outputs = virtualMas.get_outputs();
    OpenMagnetics::Inputs devirtualizedInputs;
    devirtualizedInputs.set_design_requirements(inputs.get_design_requirements());
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); operatingPointIndex++)
    {
        OperatingPoint devirtualizedOperatingPoint;
        devirtualizedOperatingPoint.set_conditions(inputs.get_operating_points()[operatingPointIndex].get_conditions());
        for (size_t windingIndex = 0; windingIndex < originalMagnetic.get_coil().get_functional_description().size(); windingIndex++)
        {
            bool found = false;
            for (size_t virtualWindingIndex = 0; virtualWindingIndex < virtualMagnetic.get_coil().get_functional_description().size(); virtualWindingIndex++)
            {
                auto virtualWindingName = virtualMagnetic.get_coil().get_functional_description()[virtualWindingIndex].get_name();
                auto windingIndexes = devirtualizationData[virtualWindingName].get_mutable_winding_indexes();
                if (std::find(windingIndexes.begin(), windingIndexes.end(), windingIndex) != windingIndexes.end()) {
                    found = true;
                    if (devirtualizationData[virtualWindingName].get_is_parallel()) {
                        double proportion = double(originalMagnetic.get_coil().get_functional_description()[windingIndex].get_number_turns()) / devirtualizationData[virtualWindingName].get_number_turns();
                        auto devirtualizedExcitation = OpenMagnetics::Inputs::get_excitation_with_proportional_voltage(inputs.get_operating_points()[operatingPointIndex].get_excitations_per_winding()[virtualWindingIndex], proportion);
                        devirtualizedOperatingPoint.get_mutable_excitations_per_winding().push_back(devirtualizedExcitation);
                    }
                    else {
                        double proportion = double(originalMagnetic.get_coil().get_functional_description()[windingIndex].get_number_turns()) / devirtualizationData[virtualWindingName].get_number_turns();
                        auto devirtualizedExcitation = OpenMagnetics::Inputs::get_excitation_with_proportional_current(inputs.get_operating_points()[operatingPointIndex].get_excitations_per_winding()[virtualWindingIndex], proportion);
                        devirtualizedOperatingPoint.get_mutable_excitations_per_winding().push_back(devirtualizedExcitation);
                    }
                    break;
                }
            }
            if (!found)
            {
                MAS::OperatingPointExcitation emptyExcitation;
                auto magnetizingInductance = OpenMagnetics::MagnetizingInductance();
                auto magnetizingInductanceOutput = magnetizingInductance.calculate_inductance_from_number_turns_and_gapping(originalMagnetic.get_core(), originalMagnetic.get_coil(), &inputs.get_mutable_operating_points()[operatingPointIndex]);
                auto magnetizingInductanceValue = OpenMagnetics::resolve_dimensional_values(magnetizingInductanceOutput.get_magnetizing_inductance());
                emptyExcitation.set_frequency(devirtualizedOperatingPoint.get_excitations_per_winding()[0].get_frequency());
                MAS::SignalDescriptor current;
                auto currentWaveform = OpenMagnetics::Inputs::create_waveform(WaveformLabel::SINUSOIDAL, 1, devirtualizedOperatingPoint.get_excitations_per_winding()[0].get_frequency());
                current.set_waveform(currentWaveform);
                emptyExcitation.set_current(current);

                MAS::SignalDescriptor voltage;
                auto voltageWaveform = OpenMagnetics::Inputs::create_waveform(WaveformLabel::SINUSOIDAL, 1, devirtualizedOperatingPoint.get_excitations_per_winding()[0].get_frequency());
                voltage.set_waveform(currentWaveform);
                emptyExcitation.set_voltage(voltage);

                auto magnetizingCurrent = OpenMagnetics::Inputs::calculate_magnetizing_current(emptyExcitation, magnetizingInductanceValue, false, 0.0);
                emptyExcitation.set_magnetizing_current(magnetizingCurrent);

                devirtualizedOperatingPoint.get_mutable_excitations_per_winding().push_back(emptyExcitation);
            }

        }
        devirtualizedInputs.get_mutable_operating_points().push_back(devirtualizedOperatingPoint);
    }

    OpenMagnetics::Mas devirtualizedMas;
    devirtualizedMas.set_inputs(devirtualizedInputs);
    devirtualizedMas.set_magnetic(originalMagnetic);
    devirtualizedMas.set_outputs(outputs);

    return devirtualizedMas;
}

std::vector<MagneticCombinator::WindingCombination> MagneticCombinator::calculate_coil_combinations(OpenMagnetics::Coil coil, int numberDesiredWindings) {
    std::vector<MagneticCombinator::WindingCombination> combinations;
    if (coil.get_functional_description().size() < size_t(numberDesiredWindings))
    {
        return combinations;
    }

    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); windingIndex++)
    {
        MagneticCombinator::WindingCombination combination;
        // combination.set_isolation_side(coil.get_functional_description()[windingIndex].get_isolation_side() == IsolationSide::PRIMARY? IsolationSide::PRIMARY : IsolationSide::SECONDARY);
        combination.set_isolation_side(coil.get_functional_description()[windingIndex].get_isolation_side());
        combination.get_mutable_winding_indexes().push_back(windingIndex);
        combination.set_connections(coil.get_functional_description()[windingIndex].get_connections());
        combination.set_number_turns(coil.get_functional_description()[windingIndex].get_number_turns());
        combination.set_number_parallels(coil.get_functional_description()[windingIndex].get_number_parallels());
        combination.set_name(coil.get_functional_description()[windingIndex].get_name());
        combinations.push_back(combination);
    }

    for (size_t windingIndex = 0; windingIndex < std::min(maximumNumberCombinationLevels, coil.get_functional_description().size() - numberDesiredWindings); windingIndex++)
    {
        auto aux = calculate_combinations(combinations, coil.get_functional_description().size() - numberDesiredWindings + 1);
        for (auto elem : aux) {
            combinations.push_back(elem);
        }
    }

    return combinations;
}

std::vector<MagneticCombinator::WindingCombination> MagneticCombinator::calculate_combinations(std::vector<MagneticCombinator::WindingCombination> sourceCombinations, int maximumWindings) {
    std::vector<MagneticCombinator::WindingCombination> combinations;
    for (size_t firstCombinationIndex = 0; firstCombinationIndex < sourceCombinations.size(); firstCombinationIndex++)
    {
        for (size_t secondCombinationIndex = firstCombinationIndex + 1; secondCombinationIndex < sourceCombinations.size(); secondCombinationIndex++)
        {
            if (repeated_windings(sourceCombinations[firstCombinationIndex].get_winding_indexes(), sourceCombinations[secondCombinationIndex].get_winding_indexes()))
                continue;

            if (sourceCombinations[firstCombinationIndex].get_isolation_side() != sourceCombinations[secondCombinationIndex].get_isolation_side()) {
                continue;
            }

            // Auxiliaries are marked as tertiary
            if (sourceCombinations[firstCombinationIndex].get_isolation_side() == IsolationSide::TERTIARY || sourceCombinations[secondCombinationIndex].get_isolation_side() == IsolationSide::TERTIARY) {
                continue;
            }
            // Shielding are marked as quaternary
            if (sourceCombinations[firstCombinationIndex].get_isolation_side() == IsolationSide::QUATERNARY || sourceCombinations[secondCombinationIndex].get_isolation_side() == IsolationSide::QUATERNARY) {
                continue;
            }

            auto compulsorySeriesConnections = get_compulsory_series_connections(firstCombinationIndex, sourceCombinations);
            auto isSecondCombinationIndexCompulsory = false;
            for (auto compulsorySeriesConnection : compulsorySeriesConnections) {
                if (sourceCombinations[secondCombinationIndex].get_name() == compulsorySeriesConnection) {
                    isSecondCombinationIndexCompulsory = true;
                }
            }

            if (compulsorySeriesConnections.size() > 0 && !isSecondCombinationIndexCompulsory) {
                // This means that the winding is tapped with another one, and that other one is the only possible combination.
                continue;
            }

            if ((sourceCombinations[firstCombinationIndex].get_winding_indexes().size() + sourceCombinations[secondCombinationIndex].get_winding_indexes().size()) <= size_t(maximumWindings)) {
                MagneticCombinator::WindingCombination combination;
                for (auto elem : sourceCombinations[firstCombinationIndex].get_winding_indexes()) {
                    combination.get_mutable_winding_indexes().push_back(elem);
                }
                for (auto elem : sourceCombinations[secondCombinationIndex].get_winding_indexes()) {
                    combination.get_mutable_winding_indexes().push_back(elem);
                }

                if (sourceCombinations[firstCombinationIndex].get_connections() && sourceCombinations[secondCombinationIndex].get_connections()) {
                    auto firstConnections = sourceCombinations[firstCombinationIndex].get_connections().value();
                    auto secondConnections = sourceCombinations[secondCombinationIndex].get_connections().value();

                    for (size_t connectionIndex = 0; connectionIndex < secondConnections.size(); connectionIndex++) {
                        firstConnections.push_back(secondConnections[connectionIndex]);
                    }
                    combination.set_connections(firstConnections);
                }
                combination.set_name(sourceCombinations[firstCombinationIndex].get_name() + "-" + sourceCombinations[secondCombinationIndex].get_name());
                combination.set_isolation_side(sourceCombinations[firstCombinationIndex].get_isolation_side());
                combination.set_number_turns(sourceCombinations[firstCombinationIndex].get_number_turns() + sourceCombinations[secondCombinationIndex].get_number_turns());
                combination.set_number_parallels(std::max(sourceCombinations[firstCombinationIndex].get_number_parallels(), sourceCombinations[secondCombinationIndex].get_number_parallels()));
                combination.set_is_parallel(false);
                combination.get_mutable_winding_combinations().push_back(sourceCombinations[firstCombinationIndex]);
                combination.get_mutable_winding_combinations().push_back(sourceCombinations[secondCombinationIndex]);
                combinations.push_back(combination);
            }
            if (sourceCombinations[firstCombinationIndex].get_number_turns() == sourceCombinations[secondCombinationIndex].get_number_turns() &&
                !sourceCombinations[secondCombinationIndex].get_is_parallel() && !sourceCombinations[secondCombinationIndex].get_is_parallel()) {


                MagneticCombinator::WindingCombination combination;
                for (auto elem : sourceCombinations[firstCombinationIndex].get_winding_indexes()) {
                    combination.get_mutable_winding_indexes().push_back(elem);
                }
                for (auto elem : sourceCombinations[secondCombinationIndex].get_winding_indexes()) {
                    combination.get_mutable_winding_indexes().push_back(elem);
                }

                if (sourceCombinations[firstCombinationIndex].get_connections() && sourceCombinations[secondCombinationIndex].get_connections()) {
                    auto firstConnections = sourceCombinations[firstCombinationIndex].get_connections().value();
                    auto secondConnections = sourceCombinations[secondCombinationIndex].get_connections().value();

                    for (size_t connectionIndex = 0; connectionIndex < secondConnections.size(); connectionIndex++) {
                        firstConnections.push_back(secondConnections[connectionIndex]);
                    }
                    combination.set_connections(firstConnections);
                }
                combination.set_name(sourceCombinations[firstCombinationIndex].get_name() + "-" + sourceCombinations[secondCombinationIndex].get_name());
                combination.set_isolation_side(sourceCombinations[firstCombinationIndex].get_isolation_side());
                combination.set_number_turns(sourceCombinations[firstCombinationIndex].get_number_turns());
                combination.set_number_parallels(sourceCombinations[firstCombinationIndex].get_number_parallels() + sourceCombinations[secondCombinationIndex].get_number_parallels());
                combination.set_is_parallel(true);
                combination.get_mutable_winding_combinations().push_back(sourceCombinations[firstCombinationIndex]);
                combination.get_mutable_winding_combinations().push_back(sourceCombinations[secondCombinationIndex]);
                combinations.push_back(combination);
            }
        }
    }

    return combinations;
}

} // namespace OpenMagnetics
