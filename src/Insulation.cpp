#include "Utils.h"
#include "Insulation.h"
#include "Settings.h"
#include <magic_enum.hpp>
#include "WindingSkinEffectLosses.h"
#include <cfloat>


namespace OpenMagnetics {

double linear_table_interpolation(std::vector<std::pair<double, double>> table, double x){
    if (x > table.back().first) {
        double slope = (table[table.size() - 1].second - table[table.size() - 2].second) / (table[table.size() - 1].first - table[table.size() - 2].first);
        return (x - table.back().first) * slope + table.back().second;
    }
    if (x < table.front().first) {
        double slope = (table[1].second - table[0].second) / (table[1].first - table[0].first);
        return (table[0].first - x) * slope + table.front().second;
    }


    for (size_t i = 0; i < table.size() - 1; ++i)
    {
        if (table[i].first <= x && x <= table[i + 1].first) {
            double proportion = (x - table[i].first) / (table[i + 1].first - table[i].first);
            return std::lerp(table[i].second, table[i + 1].second, proportion);
        }
    }

    return DBL_MAX;
}

InsulationCoordinationOutput InsulationCoordinator::calculate_insulation_coordination(InputsWrapper& inputs) {
    InsulationCoordinationOutput insulationCoordinationOutput;
    insulationCoordinationOutput.set_clearance(calculate_clearance(inputs));
    insulationCoordinationOutput.set_creepage_distance(calculate_creepage_distance(inputs, true));
    insulationCoordinationOutput.set_withstand_voltage(calculate_withstand_voltage(inputs));
    insulationCoordinationOutput.set_distance_through_insulation(calculate_distance_through_insulation(inputs));

    return insulationCoordinationOutput;
}


bool InsulationCoordinator::can_fully_insulated_wire_be_used(InputsWrapper& inputs) {
    auto standards = inputs.get_standards();
    for (auto standard : standards) {
        if (standard == InsulationStandards::IEC_603351 || standard == InsulationStandards::IEC_606641) {
            return false;
        }
    }
    return true;
}

size_t times_withstand_voltage_is_covered_by_wires(WireWrapper leftWire, WireWrapper rightWire, double withstandVoltage, bool canFullyInsulatedWireBeUsed) {
    auto leftCoatingMaybe = leftWire.resolve_coating();
    auto rightCoatingMaybe = rightWire.resolve_coating();

    size_t times = 0;
    if (leftCoatingMaybe) {
        auto coating = leftCoatingMaybe.value();

        if (leftWire.get_type() == WireType::LITZ && !coating.get_breakdown_voltage()) {
            auto strand = leftWire.resolve_strand();
            coating = WireWrapper::resolve_coating(strand).value();
        }

        if (!coating.get_breakdown_voltage()) {
            throw std::runtime_error("Wire " + leftWire.get_name().value() + " is missing breakdown voltage");
        }
        if (coating.get_breakdown_voltage().value() > withstandVoltage) {
            if (coating.get_number_layers()) {
                times += coating.get_number_layers().value();
            }
            else if (coating.get_grade() && canFullyInsulatedWireBeUsed) {
                if (coating.get_grade().value() > 3) 
                times += 3;
            }
            else {
                times++;
            }
        }
    }
    if (rightCoatingMaybe) {
        auto coating = rightCoatingMaybe.value();

        if (rightWire.get_type() == WireType::LITZ && !coating.get_breakdown_voltage()) {
            auto strand = rightWire.resolve_strand();
            coating = WireWrapper::resolve_coating(strand).value();
        }

        if (!coating.get_breakdown_voltage()) {
            throw std::runtime_error("Wire " + rightWire.get_name().value() + " is missing breakdown voltage");
        }

        if (coating.get_breakdown_voltage().value() > withstandVoltage) {
            if (coating.get_number_layers()) {
                times += coating.get_number_layers().value();
            }
            else if (coating.get_grade() && canFullyInsulatedWireBeUsed) {
                if (coating.get_grade().value() > 3) 
                times += 3;
            }
            else {
                times++;
            }
        }
    }

    return times;
}

double insulation_distance_provided_by_wires(WireWrapper leftWire, WireWrapper rightWire, double withstandVoltage, bool canFullyInsulatedWireBeUsed) {
    double distance = 0;
    {
        auto coating = leftWire.resolve_coating().value();

        if (leftWire.get_type() == WireType::LITZ && !coating.get_breakdown_voltage()) {
            auto strand = leftWire.resolve_strand();
            coating = WireWrapper::resolve_coating(strand).value();
        }

        if (!coating.get_breakdown_voltage()) {
            throw std::runtime_error("Wire " + leftWire.get_name().value() + " is missing breakdown voltage");
        }

        if (coating.get_breakdown_voltage().value() > withstandVoltage) {
            if (coating.get_number_layers() || (coating.get_grade() && canFullyInsulatedWireBeUsed)) {
                if (coating.get_thickness()) {
                    distance += resolve_dimensional_values(coating.get_thickness().value());
                }
                else {
                    distance += std::min(leftWire.get_maximum_outer_width() - leftWire.get_maximum_conducting_width(), leftWire.get_maximum_outer_height() - leftWire.get_maximum_conducting_height()) / 2;
                }
            }
        }
    }
    {
        auto coating = rightWire.resolve_coating().value();

        if (rightWire.get_type() == WireType::LITZ && !coating.get_breakdown_voltage()) {
            auto strand = rightWire.resolve_strand();
            coating = WireWrapper::resolve_coating(strand).value();
        }

        if (!coating.get_breakdown_voltage()) {
            throw std::runtime_error("Wire " + rightWire.get_name().value() + " is missing breakdown voltage");
        }

        if (coating.get_breakdown_voltage().value() > withstandVoltage) {
            if (coating.get_number_layers() || (coating.get_grade() && canFullyInsulatedWireBeUsed)) {
                if (coating.get_thickness()) {
                    distance += resolve_dimensional_values(coating.get_thickness().value());
                }
                else {
                    distance += std::min(rightWire.get_maximum_outer_width() - rightWire.get_maximum_conducting_width(), rightWire.get_maximum_outer_height() - rightWire.get_maximum_conducting_height()) / 2;
                }
            }
        }
    }

    return distance;
}

std::optional<CoilSectionInterface> InsulationCoordinator::calculate_coil_section_interface_layers(InputsWrapper& inputs, WireWrapper leftWire, WireWrapper rightWire, InsulationMaterialWrapper insulationMaterial) {
    auto settings = Settings::GetInstance();
    bool allowMarginTape = settings->get_coil_allow_margin_tape();
    bool allowInsulatedWire = settings->get_coil_allow_insulated_wire();

    CoilSectionInterface coilSectionInterface;
    size_t numberInsulationLayers = 0;
    double tapeThickness = 0;
    double tapeDielectricStrength = 0;
    auto insulationType = inputs.get_insulation_type();
    auto canFullyInsulatedWireBeUsed = can_fully_insulated_wire_be_used(inputs);

    {
        auto aux = insulationMaterial.get_thicker_tape();
        tapeThickness = aux.first;
        tapeDielectricStrength = aux.second;
    }


    coilSectionInterface.set_layer_purpose(CoilSectionInterface::LayerPurpose::INSULATING);
    auto insulationCoordination = calculate_insulation_coordination(inputs);
    auto withstandVoltage = insulationCoordination.get_withstand_voltage();
    auto clearanceAndCreepageDistance = insulationCoordination.get_creepage_distance();
    double minimumDistanceThroughInsulation = insulationCoordination.get_distance_through_insulation();

    if (!allowMarginTape && !allowInsulatedWire) {
        throw std::runtime_error("One of the options {allowMarginTape, allowInsulatedWire} must be allowed");
    }

    switch (insulationType) {
        case InsulationType::FUNCTIONAL: {
                auto timesVoltageIsCovered = times_withstand_voltage_is_covered_by_wires(leftWire, rightWire, withstandVoltage, canFullyInsulatedWireBeUsed);
                if (timesVoltageIsCovered > 0 && allowInsulatedWire) {
                    if (clearanceAndCreepageDistance > 0 && timesVoltageIsCovered >= 3) {
                        numberInsulationLayers = 1;
                        coilSectionInterface.set_layer_purpose(CoilSectionInterface::LayerPurpose::MECHANICAL);
                        clearanceAndCreepageDistance = 0;
                    }
                    else if (clearanceAndCreepageDistance > 0 && timesVoltageIsCovered == 2) {
                        numberInsulationLayers = std::max(1.0, ceil(roundFloat(withstandVoltage / tapeDielectricStrength / tapeThickness, 1)));
                        coilSectionInterface.set_layer_purpose(CoilSectionInterface::LayerPurpose::MECHANICAL);
                        clearanceAndCreepageDistance = 0;
                    }
                    else {
                        numberInsulationLayers = 1;
                    }
                }
                else {
                    numberInsulationLayers = std::max(1.0, ceil(roundFloat(withstandVoltage / tapeDielectricStrength / tapeThickness, 1)));
                }
                break;
            }
        case InsulationType::BASIC: {
                auto timesVoltageIsCovered = times_withstand_voltage_is_covered_by_wires(leftWire, rightWire, withstandVoltage, canFullyInsulatedWireBeUsed);
                if (timesVoltageIsCovered > 0 && allowInsulatedWire) {
                    if (clearanceAndCreepageDistance > 0 && timesVoltageIsCovered >= 3) {
                        numberInsulationLayers = 1;
                        coilSectionInterface.set_layer_purpose(CoilSectionInterface::LayerPurpose::MECHANICAL);
                        clearanceAndCreepageDistance = 0;
                    }
                    else if (clearanceAndCreepageDistance > 0 && timesVoltageIsCovered == 2) {
                        numberInsulationLayers = std::max(1.0, ceil(roundFloat(withstandVoltage / tapeDielectricStrength / tapeThickness, 1)));
                        clearanceAndCreepageDistance = 0;
                    }
                    else {
                        numberInsulationLayers = 1;
                    }
                }
                else {
                    numberInsulationLayers = std::max(1.0, ceil(roundFloat(withstandVoltage / tapeDielectricStrength / tapeThickness, 1)));
                }

                break;
            }
        case InsulationType::SUPPLEMENTARY: {
                auto timesVoltageIsCovered = times_withstand_voltage_is_covered_by_wires(leftWire, rightWire, withstandVoltage, canFullyInsulatedWireBeUsed);
                if (timesVoltageIsCovered > 0 && allowInsulatedWire) {
                    if (clearanceAndCreepageDistance > 0 && timesVoltageIsCovered >= 3) {
                        numberInsulationLayers = 1;
                        coilSectionInterface.set_layer_purpose(CoilSectionInterface::LayerPurpose::MECHANICAL);
                        clearanceAndCreepageDistance = 0;
                    }
                    else if (clearanceAndCreepageDistance > 0 && timesVoltageIsCovered == 2) {
                        numberInsulationLayers = std::max(1.0, ceil(roundFloat(withstandVoltage / tapeDielectricStrength / tapeThickness, 1)));
                        clearanceAndCreepageDistance = 0;
                    }
                    else {
                        numberInsulationLayers = 1;
                    }
                }
                else {
                    numberInsulationLayers = std::max(1.0, ceil(roundFloat(withstandVoltage / tapeDielectricStrength / tapeThickness, 1)));
                }
                break;
            }
        case InsulationType::DOUBLE: {
                size_t numberInsulationLayersTogether;
                size_t numberInsulationLayersSeparated;

                auto timesWithstandVoltageIsCoveredByWires = times_withstand_voltage_is_covered_by_wires(leftWire, rightWire, withstandVoltage, canFullyInsulatedWireBeUsed);
                if (timesWithstandVoltageIsCoveredByWires >= 3 && allowInsulatedWire) {
                    numberInsulationLayersTogether = 1;
                    coilSectionInterface.set_layer_purpose(CoilSectionInterface::LayerPurpose::MECHANICAL);
                    clearanceAndCreepageDistance = 0;
                }
                else if (clearanceAndCreepageDistance > 0 && timesWithstandVoltageIsCoveredByWires == 2 && allowInsulatedWire) {
                    numberInsulationLayersTogether = std::max(1.0, ceil(roundFloat(withstandVoltage / tapeDielectricStrength / tapeThickness, 1)));
                    clearanceAndCreepageDistance = 0;
                }
                else {
                    numberInsulationLayersTogether = std::max(1.0, ceil(roundFloat(withstandVoltage / tapeDielectricStrength / tapeThickness, 1)));
                }

                auto insulation = inputs.get_mutable_design_requirements().get_insulation().value();
                insulation.set_insulation_type(InsulationType::BASIC);
                inputs.get_mutable_design_requirements().set_insulation(insulation);
                auto withstandVoltageBasic = calculate_withstand_voltage(inputs);
                insulation.set_insulation_type(InsulationType::SUPPLEMENTARY);
                inputs.get_mutable_design_requirements().set_insulation(insulation);
                auto withstandVoltageSupplementary = calculate_withstand_voltage(inputs);
                double withstandVoltageSeparated = std::max(withstandVoltageBasic, withstandVoltageSupplementary);
                timesWithstandVoltageIsCoveredByWires = times_withstand_voltage_is_covered_by_wires(leftWire, rightWire, withstandVoltageSeparated, canFullyInsulatedWireBeUsed);
                if (timesWithstandVoltageIsCoveredByWires >= 3 && allowInsulatedWire) {
                    numberInsulationLayersSeparated = 1;
                    coilSectionInterface.set_layer_purpose(CoilSectionInterface::LayerPurpose::MECHANICAL);
                    clearanceAndCreepageDistance = 0;
                }
                else if (clearanceAndCreepageDistance > 0 && timesWithstandVoltageIsCoveredByWires == 2 && allowInsulatedWire) {
                    numberInsulationLayersSeparated = std::max(1.0, ceil(roundFloat(withstandVoltage / tapeDielectricStrength / tapeThickness, 1)));
                    clearanceAndCreepageDistance = 0;
                }
                else {
                    numberInsulationLayersSeparated = std::max(1.0, ceil(roundFloat(withstandVoltage / tapeDielectricStrength / tapeThickness, 1)));
                }

                numberInsulationLayers = std::min(numberInsulationLayersTogether, numberInsulationLayersSeparated);

                break;
            }
        case InsulationType::REINFORCED: {
                auto timesWithstandVoltageIsCoveredByWires = times_withstand_voltage_is_covered_by_wires(leftWire, rightWire, withstandVoltage, canFullyInsulatedWireBeUsed);
                if (timesWithstandVoltageIsCoveredByWires >= 3 && allowInsulatedWire) {
                    numberInsulationLayers = 1;
                    coilSectionInterface.set_layer_purpose(CoilSectionInterface::LayerPurpose::MECHANICAL);
                    clearanceAndCreepageDistance = 0;
                }
                else {
                    numberInsulationLayers = std::max(1.0, ceil(roundFloat(withstandVoltage / tapeDielectricStrength / tapeThickness, 1)));
                }
                break;
            }
    }
    if (tapeThickness * numberInsulationLayers <= minimumDistanceThroughInsulation) {
        numberInsulationLayers = ceil(roundFloat(minimumDistanceThroughInsulation / tapeThickness, 1));
    }

    if (clearanceAndCreepageDistance > 0 && !allowMarginTape) {
        return std::nullopt;
    }

    coilSectionInterface.set_number_layers_insulation(numberInsulationLayers);
    coilSectionInterface.set_solid_insulation_thickness(tapeThickness * numberInsulationLayers);
    coilSectionInterface.set_total_margin_tape_distance(clearanceAndCreepageDistance);
    return coilSectionInterface;
}

double InsulationCoordinator::calculate_withstand_voltage(InputsWrapper& inputs) {
    if (!inputs.get_design_requirements().get_insulation()) {
        return 0;
    }
    double solidInsulation = 0;
    for (auto standard : inputs.get_standards()) {
        switch (standard) {
            case InsulationStandards::IEC_603351: {
                solidInsulation = std::max(solidInsulation, _insulationIEC60335Model->calculate_withstand_voltage(inputs));
                break;
            }
            case InsulationStandards::IEC_606641: {
                solidInsulation = std::max(solidInsulation, _insulationIEC60664Model->calculate_withstand_voltage(inputs));
                break;
            }
            case InsulationStandards::IEC_615581: {
                solidInsulation = std::max(solidInsulation, _insulationIEC61558Model->calculate_withstand_voltage(inputs));
                break;
            }
            case InsulationStandards::IEC_623681: {
                solidInsulation = std::max(solidInsulation, _insulationIEC62368Model->calculate_withstand_voltage(inputs));
                break;
            }
        }
    }
    return solidInsulation;
}

double InsulationCoordinator::calculate_clearance(InputsWrapper& inputs) {
    if (!inputs.get_design_requirements().get_insulation()) {
        return 0;
    }

    double clearance = 0;
    for (auto standard : inputs.get_standards()) {
        switch (standard) {
            case InsulationStandards::IEC_603351: {
                clearance = std::max(clearance, _insulationIEC60335Model->calculate_clearance(inputs));
                break;
            }
            case InsulationStandards::IEC_606641: {
                clearance = std::max(clearance, _insulationIEC60664Model->calculate_clearance(inputs));
                break;
            }
            case InsulationStandards::IEC_615581: {
                clearance = std::max(clearance, _insulationIEC61558Model->calculate_clearance(inputs));
                break;
            }
            case InsulationStandards::IEC_623681: {
                clearance = std::max(clearance, _insulationIEC62368Model->calculate_clearance(inputs));
                break;
            }
        }
    }
    return clearance;
}

double InsulationCoordinator::calculate_creepage_distance(InputsWrapper& inputs, bool includeClearance) {
    if (!inputs.get_design_requirements().get_insulation()) {
        return 0;
    }
    double creepageDistance = 0;
    for (auto standard : inputs.get_standards()) {
        switch (standard) {
            case InsulationStandards::IEC_603351: {
                creepageDistance = std::max(creepageDistance, _insulationIEC60335Model->calculate_creepage_distance(inputs, includeClearance));
                break;
            }
            case InsulationStandards::IEC_606641: {
                creepageDistance = std::max(creepageDistance, _insulationIEC60664Model->calculate_creepage_distance(inputs, includeClearance));
                break;
            }
            case InsulationStandards::IEC_615581: {
                creepageDistance = std::max(creepageDistance, _insulationIEC61558Model->calculate_creepage_distance(inputs, includeClearance));
                break;
            }
            case InsulationStandards::IEC_623681: {
                creepageDistance = std::max(creepageDistance, _insulationIEC62368Model->calculate_creepage_distance(inputs, includeClearance));
                break;
            }
        }
    }
    return creepageDistance;
}

double InsulationCoordinator::calculate_distance_through_insulation(InputsWrapper& inputs) {
    if (!inputs.get_design_requirements().get_insulation()) {
        return 0;
    }
    double dti = 0;
    for (auto standard : inputs.get_standards()) {
        switch (standard) {
            case InsulationStandards::IEC_603351: {
                dti = std::max(dti, _insulationIEC60335Model->calculate_distance_through_insulation(inputs));
                break;
            }
            case InsulationStandards::IEC_606641: {
                dti = std::max(dti, _insulationIEC60664Model->calculate_distance_through_insulation(inputs));
                break;
            }
            case InsulationStandards::IEC_615581: {
                dti = std::max(dti, _insulationIEC61558Model->calculate_distance_through_insulation(inputs));
                break;
            }
            case InsulationStandards::IEC_623681: {
                dti = std::max(dti, _insulationIEC62368Model->calculate_distance_through_insulation(inputs));
                break;
            }
        }
    }
    return dti;
}

bool InsulationIEC60664Model::electric_field_strenth_is_valid(double dti, double voltage) {
    if (dti == 0) {
        return false;
    }
    else if (dti < 30e-6) {
        return (voltage / dti) < 10e6;
    }
    else if (dti > 0.00075) {
        return (voltage / dti) < 2e6;
    }
    else {
        return (voltage / dti) < (0.25 / (dti * 1000) + 1.667) * 1e6;
    }
}

double InsulationIEC60664Model::calculate_distance_through_insulation_over_30kHz(double workingVoltage) {
    double dti = 0;
    while (!electric_field_strenth_is_valid(dti, workingVoltage)) {
        dti += 1e-6;
    }

    return dti;
}

double InsulationIEC60664Model::get_rated_impulse_withstand_voltage(OvervoltageCategory overvoltageCategory, double ratedVoltage, InsulationType insulationType){
    std::string overvoltageCategoryString = std::string{magic_enum::enum_name(overvoltageCategory)};
    auto aux = part1TableF1[overvoltageCategoryString];
    for (size_t voltagesIndex = 0; voltagesIndex < aux.size(); voltagesIndex++) {
        if (ratedVoltage <= aux[voltagesIndex].first) {
            if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
                if (voltagesIndex < aux.size() - 1) {
                    return aux[voltagesIndex + 1].second;
                }
                else {
                    return aux[voltagesIndex].second * 1.6;
                }
            }
            else {
                return aux[voltagesIndex].second;
            }
        }
    }
    throw std::invalid_argument("Too much voltage for IEC 60664-1: " + std::to_string(ratedVoltage));
}

double InsulationIEC60664Model::get_clearance_table_f2(PollutionDegree pollutionDegree, double ratedImpulseWithstandVoltage){
    std::string pollutionDegreeString = std::string{magic_enum::enum_name(pollutionDegree)};
    auto aux = part1TableF2["inhomogeneusField"][pollutionDegreeString];
    for (size_t voltagesIndex = 0; voltagesIndex < aux.size(); voltagesIndex++) {
        if (ratedImpulseWithstandVoltage <= aux[voltagesIndex].first) {
            return aux[voltagesIndex].second;
        }
    }
    throw std::invalid_argument("Too much voltage for IEC 60664-1: " + std::to_string(ratedImpulseWithstandVoltage));
}

double InsulationIEC60664Model::get_clearance_table_f8(double ratedImpulseWithstandVoltage){
    auto aux = part1TableF8["inhomogeneusField"];
    double ratedImpulseWithstandVoltageScaled = ratedImpulseWithstandVoltage;
    for (size_t voltagesIndex = 0; voltagesIndex < aux.size(); voltagesIndex++) {
        if (ratedImpulseWithstandVoltageScaled <= aux[voltagesIndex].first) {
            return aux[voltagesIndex].second;
        }
    }
    throw std::invalid_argument("Too much voltage for IEC 60664-1: " + std::to_string(ratedImpulseWithstandVoltageScaled));
}

std::optional<double> InsulationIEC60664Model::get_clearance_planar(double altitude, double ratedImpulseWithstandVoltage){
    std::vector<std::pair<double, double>> table;
    if (altitude <= lowerAltitudeLimit) {
        table = part5Table2["inhomogeneusField"];
    }
    else {
        table = part5Table3["inhomogeneusField"];
    }

    bool insideTable = false;
    for (auto& voltagePair : table) {
        if (ratedImpulseWithstandVoltage < voltagePair.first) {
            insideTable = true;
        }
    }

    if (insideTable) {
        return linear_table_interpolation(table, ratedImpulseWithstandVoltage);
    }

    // If we find to option in this table, we go back to IEC 60664-1
    return std::nullopt;
}

double InsulationIEC60664Model::get_rated_insulation_voltage(double mainSupplyVoltage){
    for (auto& voltagePair : part1TableF3) {
        if (mainSupplyVoltage < voltagePair.first) {
            return voltagePair.second;
        }
    }
    throw std::invalid_argument("Too much voltage for IEC 60664-1: " + std::to_string(mainSupplyVoltage));
}

double InsulationIEC60664Model::get_creepage_distance(PollutionDegree pollutionDegree, Cti cti, double voltageRms, WiringTechnology wiringType){
    std::string pollutionDegreeString = std::string{magic_enum::enum_name(pollutionDegree)};
    std::string ctiString = std::string{magic_enum::enum_name(cti)};
    std::string wiringTypeString = std::string{magic_enum::enum_name(wiringType)};

    if (!part1TableF5.contains(wiringTypeString)) 
        throw std::invalid_argument("Unknown wiring type: " + wiringTypeString);
    if (!part1TableF5[wiringTypeString].contains(pollutionDegreeString)) 
        throw std::invalid_argument("Pollution degree " + pollutionDegreeString + " is not supported for wiring " + wiringTypeString + " in IEC 60664");
    if (!part1TableF5[wiringTypeString][pollutionDegreeString].contains(ctiString)) 
        throw std::invalid_argument("CTI " + ctiString + " is not supported for pollution degree " + pollutionDegreeString + " and wiring " + wiringTypeString + " in IEC 60664");

    auto aux = part1TableF5[wiringTypeString][pollutionDegreeString][ctiString];
    for (auto& voltagePair : aux) {
        if (voltageRms < voltagePair.first) {
            return voltagePair.second;
        }
    }
    throw std::invalid_argument("Too much voltage for IEC 60664-1: " + std::to_string(voltageRms));
}

double InsulationIEC60664Model::get_creepage_distance_over_30kHz(double voltageRms, double frequency){
    double previousStandardFrequency = iec60664Part1MaximumFrequency;
    std::vector<std::pair<double, double>> previousTable;
    for (auto const& [standardFrequency, voltageList] : part4Table2)
    {   
        double topValue = DBL_MAX;
        if (frequency >= previousStandardFrequency && frequency <= standardFrequency) {
            topValue = linear_table_interpolation(voltageList, voltageRms);
            if (previousTable.empty()) {
                return topValue;
            }
            else {
                auto bottomValue = linear_table_interpolation(previousTable, voltageRms);
                double proportion = (frequency - previousStandardFrequency) / (standardFrequency - previousStandardFrequency);
                return std::lerp(bottomValue, topValue, proportion);
            }
        }

        previousStandardFrequency = standardFrequency;
        previousTable = voltageList;
    }
    throw std::invalid_argument("Too much frequency for IEC 60664-4: " + std::to_string(frequency));
}

std::optional<double> InsulationIEC60664Model::get_creepage_distance_planar(PollutionDegree pollutionDegree, Cti cti, double voltageRms){
    std::string pollutionDegreeString = std::string{magic_enum::enum_name(pollutionDegree)};
    std::string ctiString = std::string{magic_enum::enum_name(cti)};

    if (!part5Table4.contains(pollutionDegreeString)) 
        throw std::invalid_argument("Pollution degree " + pollutionDegreeString + " is not supported in IEC 60664-5");
    if (!part5Table4[pollutionDegreeString].contains(ctiString)) 
        throw std::invalid_argument("CTI " + ctiString + " is not supported for pollution degree " + pollutionDegreeString + " in IEC 60664-5");

    auto table = part5Table4[pollutionDegreeString][ctiString];

    bool insideTable = false;
    for (auto& voltagePair : table) {
        if (voltageRms < voltagePair.first) {
            insideTable = true;
        }
    }

    if (insideTable) {
        return linear_table_interpolation(table, voltageRms);
    }

    // If we find to option in this table, we go back to IEC 60664-1
    return std::nullopt;
}

double InsulationIEC60664Model::get_clearance_altitude_factor_correction(double altitude){
    return linear_table_interpolation(part1TableA2, altitude);
}

double InsulationIEC60664Model::get_clearance_over_30kHz(double ratedVoltagePeak, double frequency, double currentClearance){
    auto windingSkinEffectLossesModel = OpenMagnetics::WindingSkinEffectLosses();
    auto wireCurvature = windingSkinEffectLossesModel.calculate_skin_depth("copper", frequency, 20);
    bool isHomogeneus = wireCurvature >= currentClearance * 0.2;
    if (isHomogeneus) {
        double criticalFrequencyInMHz = 0.2 / currentClearance;
        double frequencyInMHz = frequency / 1e6;
        double minimumFrequencyInMHz = 3;
        if (frequencyInMHz < criticalFrequencyInMHz) {
            return currentClearance;
        }
        else if (frequencyInMHz > minimumFrequencyInMHz) {
            return currentClearance * 1.25;
        }
        else {
            double factor = 1 + (frequencyInMHz - criticalFrequencyInMHz) / (minimumFrequencyInMHz - criticalFrequencyInMHz) * 0.25;
            return currentClearance * factor;
        }
    }
    else {
        return linear_table_interpolation(part4Table1, ratedVoltagePeak);
    }
    return currentClearance;
}

double InsulationIEC60664Model::calculate_distance_through_insulation(InputsWrapper& inputs) {
    double maximumVoltageRms = inputs.get_maximum_voltage_rms();
    double maximumFrequency = inputs.get_maximum_frequency();
    double distanceThroughInsulation = 0;
    if (maximumFrequency > iec60664Part1MaximumFrequency) {
        distanceThroughInsulation = std::max(distanceThroughInsulation, calculate_distance_through_insulation_over_30kHz(maximumVoltageRms));
    }

    return ceilFloat(distanceThroughInsulation, 5);
}

double InsulationIEC60664Model::calculate_withstand_voltage(InputsWrapper& inputs) {
    double maximumVoltageRms = inputs.get_maximum_voltage_rms();
    auto overvoltageCategory = inputs.get_overvoltage_category();
    auto insulationType = inputs.get_insulation_type();

    double voltageDueToTransientOvervoltages = get_rated_impulse_withstand_voltage(overvoltageCategory, maximumVoltageRms, insulationType);
    double voltageDueToTemporaryWithstandOvervoltages = maximumVoltageRms + 1200;
    if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
        voltageDueToTemporaryWithstandOvervoltages *= 2;
    }
    double F1 = 1.2;  // defined in section 6.4.6.1
    double F3 = 1.25;  // defined in section 6.4.6.1
    double F4 = 1.1;  // defined in section 6.4.6.1
    double voltageDueToRecurringPeakVoltages = F1 * F4 * sqrt(2) * maximumVoltageRms;
    if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
        voltageDueToRecurringPeakVoltages *= F3;
    }
    double voltageDueToSteadyStateVoltages = inputs.get_maximum_voltage_peak();

    return std::max(voltageDueToTransientOvervoltages, std::max(voltageDueToTemporaryWithstandOvervoltages, std::max(voltageDueToRecurringPeakVoltages, voltageDueToSteadyStateVoltages)));
}

double InsulationIEC60664Model::calculate_clearance(InputsWrapper& inputs) {
    auto wiringTechnology = inputs.get_design_requirements().get_wiring_technology();
    auto pollutionDegree = inputs.get_pollution_degree();
    double maximumFrequency = inputs.get_maximum_frequency();
    double ratedVoltage = inputs.get_maximum_voltage_rms();
    auto overvoltageCategory = inputs.get_overvoltage_category();
    auto altitude = resolve_dimensional_values(inputs.get_altitude(), DimensionalValues::MAXIMUM);
    auto insulationType = inputs.get_insulation_type();
    double ratedImpulseWithstandVoltage = get_rated_impulse_withstand_voltage(overvoltageCategory, ratedVoltage, insulationType);
    double steadyStatePeakVoltage = inputs.get_maximum_voltage_peak();

    if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
        steadyStatePeakVoltage *= 1.6;
    }

    double clearanceToWithstandTransientOvervoltages;
    auto clearanceToWithstandTransientOvervoltagesPlanar = std::optional<double>{std::nullopt};
    if (wiringTechnology == WiringTechnology::PRINTED) {
        clearanceToWithstandTransientOvervoltagesPlanar = get_clearance_planar(altitude, ratedImpulseWithstandVoltage);
        if (clearanceToWithstandTransientOvervoltagesPlanar) {
            clearanceToWithstandTransientOvervoltages = clearanceToWithstandTransientOvervoltagesPlanar.value();
        }
    }
    if (wiringTechnology == WiringTechnology::WOUND || !clearanceToWithstandTransientOvervoltagesPlanar) {
        clearanceToWithstandTransientOvervoltages = get_clearance_table_f2(pollutionDegree, ratedImpulseWithstandVoltage);
    }

    double clearanceToWithstandSteadyStatePeakVoltages;
    auto clearanceToWithstandSteadyStatePeakVoltagesPlanar = std::optional<double>{std::nullopt};
    if (wiringTechnology == WiringTechnology::PRINTED) {
        clearanceToWithstandSteadyStatePeakVoltagesPlanar = get_clearance_planar(altitude, steadyStatePeakVoltage);
        if (clearanceToWithstandSteadyStatePeakVoltagesPlanar) {
            clearanceToWithstandSteadyStatePeakVoltages = clearanceToWithstandSteadyStatePeakVoltagesPlanar.value();
        }
    }
    if (wiringTechnology == WiringTechnology::WOUND || !clearanceToWithstandSteadyStatePeakVoltagesPlanar) {

        clearanceToWithstandSteadyStatePeakVoltages = get_clearance_table_f8(steadyStatePeakVoltage);
        if (maximumFrequency > iec60664Part1MaximumFrequency) {
            clearanceToWithstandSteadyStatePeakVoltages = get_clearance_over_30kHz(steadyStatePeakVoltage, maximumFrequency, clearanceToWithstandSteadyStatePeakVoltages);
        }
    }

    double clearance = std::max(clearanceToWithstandTransientOvervoltages, clearanceToWithstandSteadyStatePeakVoltages);

    if (altitude > lowerAltitudeLimit) {
        clearance *= get_clearance_altitude_factor_correction(altitude);
    }

    return clearance;
}

double InsulationIEC60664Model::calculate_creepage_distance(InputsWrapper& inputs, bool includeClearance) {
    auto wiringTechnology = inputs.get_design_requirements().get_wiring_technology();
    auto pollutionDegree = inputs.get_pollution_degree();
    auto cti = inputs.get_cti();
    auto insulationType = inputs.get_insulation_type();
    double maximumFrequency = inputs.get_maximum_frequency();
    double maximumVoltageRms = inputs.get_maximum_voltage_rms();
    double ratedVoltage = inputs.get_maximum_voltage_rms();
    double mainSupplyVoltage = resolve_dimensional_values(inputs.get_main_supply_voltage());
    double ratedInsulationVoltage = get_rated_insulation_voltage(mainSupplyVoltage);

    double voltageRms = std::max(maximumVoltageRms, std::max(ratedVoltage, ratedInsulationVoltage));

    double creepageDistance;
    auto creepageDistancePlanar = std::optional<double>{std::nullopt};
    if (wiringTechnology == WiringTechnology::PRINTED) {
        creepageDistancePlanar = get_creepage_distance_planar(pollutionDegree, cti, voltageRms);
        if (creepageDistancePlanar) {
            creepageDistance = creepageDistancePlanar.value();
        }
    }

    if (wiringTechnology == WiringTechnology::WOUND || !creepageDistancePlanar) {
        creepageDistance = get_creepage_distance(pollutionDegree, cti, voltageRms);
        if (maximumFrequency > iec60664Part1MaximumFrequency) {
            double maximumVoltagePeak = inputs.get_maximum_voltage_peak();
            double creepageDistanceOver30kHz = get_creepage_distance_over_30kHz(maximumVoltagePeak, maximumFrequency);
            switch (pollutionDegree) {
                case PollutionDegree::P1:
                    break;
                case PollutionDegree::P2:
                    creepageDistanceOver30kHz *= 1.2;  // According to table 2 of IEC 60664-4
                    break;
                case PollutionDegree::P3:
                    creepageDistanceOver30kHz *= 1.4;  // According to table 2 of IEC 60664-4
                    break;
            }
            creepageDistance = std::max(creepageDistance, creepageDistanceOver30kHz);
        }
    }

    if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
        creepageDistance *= 2;
    }


    bool allowSmallerCreepageDistanceThanClearance = false;
    if (pollutionDegree == PollutionDegree::P1 || pollutionDegree == PollutionDegree::P2) {
        auto overvoltageCategory = inputs.get_overvoltage_category();
        auto ratedImpulseWithstandVoltage = get_rated_impulse_withstand_voltage(overvoltageCategory, ratedVoltage, insulationType);
        auto clearanceToWithstandTransientOvervoltages = get_clearance_table_f2(pollutionDegree, ratedImpulseWithstandVoltage);
        if (clearanceToWithstandTransientOvervoltages < creepageDistance) {
            allowSmallerCreepageDistanceThanClearance = true;
        }
    }

    if (!allowSmallerCreepageDistanceThanClearance && includeClearance) {
        creepageDistance = std::max(creepageDistance, calculate_clearance(inputs));
    }

    return roundFloat(creepageDistance, 5);
}

double InsulationIEC62368Model::get_working_voltage(InputsWrapper& inputs) {
    return inputs.get_maximum_voltage_peak();
}

double InsulationIEC62368Model::get_working_voltage_rms(InputsWrapper& inputs) {
    return inputs.get_maximum_voltage_rms();
}

double InsulationIEC62368Model::get_required_withstand_voltage(InputsWrapper& inputs) {
    return get_working_voltage(inputs);
}

double InsulationIEC62368Model::get_voltage_due_to_temporary_overvoltages_procedure_1(double supplyVoltage){
    double voltageDueToTemporaryWithstandOvervoltages = supplyVoltage + 1200;
    if (supplyVoltage <= 250) {
        return std::max(2000.0, voltageDueToTemporaryWithstandOvervoltages);
    }
    else {
        return std::max(2500.0, voltageDueToTemporaryWithstandOvervoltages);
    }
}

double InsulationIEC62368Model::get_voltage_due_to_transient_overvoltages(double requiredWithstandVoltage, InsulationType insulationType){
    std::vector<std::pair<double, double>> table;
    if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
        table = table25["REINFORCED"];
    }
    else {
        table = table25["BASIC"];
    }
    return linear_table_interpolation(table, requiredWithstandVoltage);
}

double InsulationIEC62368Model::get_voltage_due_to_recurring_peak_voltages(double workingVoltage, InsulationType insulationType){
    std::vector<std::pair<double, double>> table;
    if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
        table = table26["REINFORCED"];
    }
    else {
        table = table26["BASIC"];
    }
    return linear_table_interpolation(table, workingVoltage);
}

double InsulationIEC62368Model::get_voltage_due_to_temporary_overvoltages(double supplyVoltageRms, InsulationType insulationType){
    std::vector<std::pair<double, double>> table;
    if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
        table = table27["REINFORCED"];
    }
    else {
        table = table27["BASIC"];
    }
    for (auto& voltagePair : table) {
        if (supplyVoltageRms < voltagePair.first) {
            return voltagePair.second;
        }
    }
    throw std::invalid_argument("Too much voltage for IEC 62368-1 in table 27: " + std::to_string(supplyVoltageRms));
}

double InsulationIEC62368Model::get_reduction_factor_per_material(std::string material, double frequency){
    double previousStandardFrequency = iec62368LowerFrequency;
    for (auto const& [standardFrequency, reductionFactor] : table22[material]) {
        if (frequency >= previousStandardFrequency && frequency <= standardFrequency) {
            return reductionFactor;
        }
        previousStandardFrequency = standardFrequency;
    }
    throw std::invalid_argument("Too much frequency for IEC 62368-1 in table 22: " + std::to_string(frequency));
}

double InsulationIEC62368Model::get_clearance_table_10(double supplyVoltagePeak, InsulationType insulationType, PollutionDegree pollutionDegree){
    std::string pollutionDegreeString = std::string{magic_enum::enum_name(pollutionDegree)};
    std::vector<std::pair<double, double>> table;
    if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
        table = table10["REINFORCED"][pollutionDegreeString];
    }
    else {
        table = table10["BASIC"][pollutionDegreeString];
    }

    double result = linear_table_interpolation(table, supplyVoltagePeak);

    if (result < 0.005) {
        return ceilFloat(result, 5);
    }
    else {
        return ceilFloat(result, 4);
    }
}

double InsulationIEC62368Model::get_clearance_table_14(double supplyVoltagePeak, InsulationType insulationType, PollutionDegree pollutionDegree){
    std::string pollutionDegreeString = std::string{magic_enum::enum_name(pollutionDegree)};
    std::vector<std::pair<double, double>> table;
    if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
        table = table14["REINFORCED"][pollutionDegreeString];
    }
    else {
        table = table14["BASIC"][pollutionDegreeString];
    }

    double result = linear_table_interpolation(table, supplyVoltagePeak);
    
    if (result < 0.005) {
        return ceilFloat(result, 5);
    }
    else {
        return ceilFloat(result, 4);
    }
}

double InsulationIEC62368Model::get_clearance_table_11(double supplyVoltagePeak, InsulationType insulationType, PollutionDegree pollutionDegree){
    std::vector<std::pair<double, double>> table;
    if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
        table = table11["REINFORCED"];
    }
    else {
        table = table11["BASIC"];
    }

    double valuePollutionDegree2 = linear_table_interpolation(table, supplyVoltagePeak);
    double result;
    if (pollutionDegree == PollutionDegree::P1) {
        result = valuePollutionDegree2 * 0.8;
    }
    else if (pollutionDegree == PollutionDegree::P3) {
        result = valuePollutionDegree2 * 1.4;
    }
    else {
        result = valuePollutionDegree2;
    }

    if (result < 0.005) {
        return ceilFloat(result, 5);
    }
    else {
        return ceilFloat(result, 4);
    }
}

double InsulationIEC62368Model::get_distance_table_G13(double workingVoltage, InsulationType insulationType){
    std::vector<std::pair<double, double>> table;
    if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
        table = tableG13["REINFORCED"];
    }
    else {
        table = tableG13["BASIC"];
    }

    double result = linear_table_interpolation(table, workingVoltage);
    return ceilFloat(result, 4);
}

double InsulationIEC62368Model::get_creepage_distance_table_17(double voltageRms, InsulationType insulationType, PollutionDegree pollutionDegree, Cti cti){
    std::string pollutionDegreeString = std::string{magic_enum::enum_name(pollutionDegree)};
    std::string ctiString = std::string{magic_enum::enum_name(cti)};
    std::vector<std::pair<double, double>> table = table17[pollutionDegreeString][ctiString];

    double valueBasic = linear_table_interpolation(table, voltageRms);

    if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
        return ceilFloat(valueBasic * 2, 4);
    }
    else {
        return ceilFloat(valueBasic, 4);
    }
}

double InsulationIEC62368Model::get_creepage_distance_table_18(double voltageRms, double frequency, PollutionDegree pollutionDegree, InsulationType insulationType){
    double previousStandardFrequency = iec62368LowerFrequency;
    double valuePollutionDegree1 = DBL_MAX;
    for (auto const& [standardFrequencyString, voltageList] : table18)
    {
        double standardFrequency = standardFrequencyString;
        if ((frequency >= previousStandardFrequency && frequency <= standardFrequency) || standardFrequency == 400000) {
            valuePollutionDegree1 = linear_table_interpolation(voltageList, voltageRms);
            break;
        }
        previousStandardFrequency = standardFrequency;
    }
    if (valuePollutionDegree1 == DBL_MAX) {
        throw std::invalid_argument("Too much frequency for IEC 62368-1 in table 18: " + std::to_string(frequency));
    }

    double result;
    if (pollutionDegree == PollutionDegree::P1) {
        result = valuePollutionDegree1;
    }
    else if (pollutionDegree == PollutionDegree::P2) {
        result = valuePollutionDegree1 * 1.2;
    }
    else {
        result = valuePollutionDegree1 * 1.4;
    }

    if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
        return ceilFloat(2 * result, 5);
    }
    else {
        return ceilFloat(result, 5);
    }
}

double InsulationIEC62368Model::get_altitude_factor(double altitude){
    double result = linear_table_interpolation(table16, altitude);
    return ceilFloat(result, 5);
}

double InsulationIEC62368Model::get_mains_transient_voltage(double supplyVoltagePeak, OvervoltageCategory overvoltageCategory){
    std::string overvoltageCategoryString = std::string{magic_enum::enum_name(overvoltageCategory)};
    auto table = table12[overvoltageCategoryString];
    for (size_t voltagesIndex = 0; voltagesIndex < table.size(); voltagesIndex++) {
        if (supplyVoltagePeak <= table[voltagesIndex].first) {
            return table[voltagesIndex].second;
        }
    }
    throw std::invalid_argument("Too much voltage for IEC 62368-1 in table 12: " + std::to_string(supplyVoltagePeak));
}

double InsulationIEC62368Model::get_es2_voltage_limit(double frequency){
    if (frequency < 1) {
        return 120;
    }
    else if (frequency < 1000) {
        return 50;
    }
    else if (frequency < 100000) {
        return 50 + 0.9 * (frequency / 1000);
    }
    else {
        return 140;
    }
}

double InsulationIEC62368Model::calculate_withstand_voltage(InputsWrapper& inputs) {
    // double maximumFrequency = inputs.get_maximum_frequency();
    double workingVoltage = get_working_voltage(inputs);
    double requiredWithstandVoltage = get_required_withstand_voltage(inputs);
    double mainSupplyVoltage = resolve_dimensional_values(inputs.get_main_supply_voltage());
    auto insulationType = inputs.get_insulation_type();

    // For when getting the breakdown voltage of a material, we need to apply this reduction factor
    // if (maximumFrequency > iec62368LowerFrequency) {
    //     auto reductionFactor = get_reduction_factor_per_material("Other thin foil materials", maximumFrequency);
    //     if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
    //         materialBreakDownVoltage = 1.2 * 2 * materialBreakDownVoltage / reductionFactor;
    //     }
    //     else {
    //         materialBreakDownVoltage = 1.2 * materialBreakDownVoltage / reductionFactor;
    //     }
    // }

    double voltageDueToTransientOvervoltages = get_voltage_due_to_transient_overvoltages(workingVoltage, insulationType);
    double voltageDueToRecurringPeakVoltages = get_voltage_due_to_recurring_peak_voltages(requiredWithstandVoltage, insulationType);
    double voltageDueToTemporaryOvervoltages = get_voltage_due_to_temporary_overvoltages(mainSupplyVoltage, insulationType);


    return std::max(voltageDueToTransientOvervoltages, std::max(voltageDueToRecurringPeakVoltages, voltageDueToTemporaryOvervoltages));
}

double InsulationIEC62368Model::calculate_distance_through_insulation(InputsWrapper& inputs) {
    double maximumFrequency = inputs.get_maximum_frequency();
    double es2VoltageLimit = get_es2_voltage_limit(maximumFrequency);
    double workingVoltageRms = get_working_voltage_rms(inputs);
    auto insulationType = inputs.get_insulation_type();

    if (workingVoltageRms <= es2VoltageLimit) {
        return 0;
    }
    else {
        if (insulationType == InsulationType::FUNCTIONAL || insulationType == InsulationType::BASIC) {
            return 0;
        }
        else {
            return 0.0004;
        }
    }
}

double InsulationIEC62368Model::calculate_clearance(InputsWrapper& inputs) {
    auto wiringTechnology = inputs.get_design_requirements().get_wiring_technology();
    auto pollutionDegree = inputs.get_pollution_degree();
    auto overvoltageCategory = inputs.get_overvoltage_category();
    double maximumFrequency = inputs.get_maximum_frequency();
    double voltagePeak = inputs.get_maximum_voltage_peak();
    double workingVoltage = get_working_voltage(inputs);
    double requiredWithstandVoltage = get_required_withstand_voltage(inputs);
    double mainSupplyVoltage = resolve_dimensional_values(inputs.get_main_supply_voltage());
    double voltageDueToTemporaryOvervoltages = get_voltage_due_to_temporary_overvoltages_procedure_1(mainSupplyVoltage);
    double voltageProcedure1 = std::max(voltagePeak, std::max(workingVoltage, voltageDueToTemporaryOvervoltages));
    auto altitude = resolve_dimensional_values(inputs.get_altitude(), DimensionalValues::MAXIMUM);
    auto insulationType = inputs.get_insulation_type();

    if (wiringTechnology == WiringTechnology::PRINTED) {
        return get_distance_table_G13(voltagePeak, insulationType);
    }

    double clearanceProcedure1;
    double clearanceProcedure2;
    if (maximumFrequency <= iec62368LowerFrequency) {
        clearanceProcedure1 = get_clearance_table_10(voltageProcedure1, insulationType, pollutionDegree);
    }
    else {
        clearanceProcedure1 = get_clearance_table_11(voltagePeak, insulationType, pollutionDegree);
    }
    double mainsTransientVoltage = get_mains_transient_voltage(mainSupplyVoltage, overvoltageCategory);
    clearanceProcedure2 = get_clearance_table_14(std::max(requiredWithstandVoltage, mainsTransientVoltage), insulationType, pollutionDegree);

    double altitudeFactor = get_altitude_factor(altitude);

    return ceilFloat(altitudeFactor * std::max(clearanceProcedure1, clearanceProcedure2), 5);
    // TODO: remove distance if FIW complies with conditions in table G.4
}

double InsulationIEC62368Model::calculate_creepage_distance(InputsWrapper& inputs, bool includeClearance) {
    auto wiringTechnology = inputs.get_design_requirements().get_wiring_technology();
    double voltagePeak = inputs.get_maximum_voltage_peak();
    double workingVoltageRms = get_working_voltage_rms(inputs);
    auto cti = inputs.get_cti();
    auto pollutionDegree = inputs.get_pollution_degree();
    double maximumFrequency = inputs.get_maximum_frequency();
    auto insulationType = inputs.get_insulation_type();

    if (wiringTechnology == WiringTechnology::PRINTED) {
        return get_distance_table_G13(voltagePeak, insulationType);
    }

    double creepageDistance = DBL_MAX;

    if (maximumFrequency <= iec62368LowerFrequency) {
        creepageDistance = get_creepage_distance_table_17(workingVoltageRms, insulationType, pollutionDegree, cti);
    }
    else {
        creepageDistance = std::max(get_creepage_distance_table_17(workingVoltageRms, insulationType, pollutionDegree, cti),
                        get_creepage_distance_table_18(workingVoltageRms, maximumFrequency, pollutionDegree, insulationType));
    }

    if (includeClearance) {
        creepageDistance = std::max(creepageDistance, calculate_clearance(inputs));
    }

    return creepageDistance;

    // TODO: remove distance if FIW complies with conditions in table G.4
}

double InsulationIEC61558Model::get_working_voltage_peak(InputsWrapper& inputs) {
    return inputs.get_maximum_voltage_peak();
}

double InsulationIEC61558Model::get_working_voltage_rms(InputsWrapper& inputs) {
    return inputs.get_maximum_voltage_rms();
}

double InsulationIEC61558Model::get_withstand_voltage_table_14(OvervoltageCategory overvoltageCategory, InsulationType insulationType, double workingVoltage){
    std::string overvoltageCategoryString = std::string{magic_enum::enum_name(overvoltageCategory)};
    std::string insulationTypeString = std::string{magic_enum::enum_name(insulationType)};
    std::vector<std::pair<double, double>> table = table14[overvoltageCategoryString][insulationTypeString];

    return linear_table_interpolation(table, workingVoltage);
}

double InsulationIEC61558Model::get_clearance_table_20(OvervoltageCategory overvoltageCategory, PollutionDegree pollutionDegree, InsulationType insulationType, double workingVoltage){
    std::string pollutionDegreeString = std::string{magic_enum::enum_name(pollutionDegree)};
    std::string overvoltageCategoryString = std::string{magic_enum::enum_name(overvoltageCategory)};
    std::string insulationTypeString = std::string{magic_enum::enum_name(insulationType)};
    std::vector<std::pair<double, double>> table = table20[overvoltageCategoryString][insulationTypeString][pollutionDegreeString];

    if (workingVoltage < iec61558MinimumWorkingVoltage || pollutionDegree == PollutionDegree::P1 || insulationType == InsulationType::FUNCTIONAL) {
        return 0;
    }

    for (size_t voltagesIndex = 0; voltagesIndex < table.size(); voltagesIndex++) {
        if (workingVoltage <= table[voltagesIndex].first) {
            return table[voltagesIndex].second;
        }
    }
    throw std::invalid_argument("Too much voltage for IEC 61558-1: " + std::to_string(workingVoltage));
}

double InsulationIEC61558Model::get_creepage_distance_table_21(Cti cti, PollutionDegree pollutionDegree, InsulationType insulationType, double workingVoltage){
    std::string pollutionDegreeString = std::string{magic_enum::enum_name(pollutionDegree)};
    std::string ctiString = std::string{magic_enum::enum_name(cti)};
    std::string insulationTypeString = std::string{magic_enum::enum_name(insulationType)};
    std::vector<std::pair<double, double>> table = table21[ctiString][insulationTypeString][pollutionDegreeString];

    if (workingVoltage < iec61558MinimumWorkingVoltage || insulationType == InsulationType::FUNCTIONAL) {
        return 0;
    }

    double creepageDistance = linear_table_interpolation(table, workingVoltage);

    return creepageDistance;
}

double InsulationIEC61558Model::get_distance_through_insulation_table_22(InsulationType insulationType, double workingVoltage, bool usingThinLayers){
    std::string insulationTypeString = std::string{magic_enum::enum_name(insulationType)};
    std::vector<std::pair<double, double>> table;
    if (workingVoltage < iec61558MinimumWorkingVoltage || insulationType == InsulationType::FUNCTIONAL || insulationType == InsulationType::BASIC) {
        return 0;
    }
    if (usingThinLayers) {
        table = table22[insulationTypeString]["thinLayers"];
    }
    else {
        table = table22[insulationTypeString]["solid"];
    }

    double dti = linear_table_interpolation(table, workingVoltage);
    return dti;
}

bool InsulationIEC61558Model::electric_field_strenth_is_valid(double dti, double voltage) {
    if (dti == 0) {
        return false;
    }
    else if (dti < 30e-6) {
        return (voltage / dti) < 10e6;
    }
    else if (dti > 0.00075) {
        return (voltage / dti) < 2e6;
    }
    else {
        return (voltage / dti) < (0.25 / (dti * 1000) + 1.667) * 1e6;
    }
}

double InsulationIEC61558Model::calculate_distance_through_insulation_over_30kHz(double workingVoltage) {
    double dti = 0;
    while (!electric_field_strenth_is_valid(dti, workingVoltage)) {
        dti += 1e-6;
    }

    return dti;
}

double InsulationIEC61558Model::calculate_clearance_over_30kHz(InsulationType insulationType, double workingVoltage) {
    std::string insulationTypeString = std::string{magic_enum::enum_name(insulationType)};

    if (workingVoltage < iec61558MinimumWorkingVoltage || insulationType == InsulationType::FUNCTIONAL) {
        return 0;
    }

    double clearance = DBL_MAX;

    {
        std::vector<std::pair<double, double>> table;
        if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
            table = table103["REINFORCED"];
        }
        else {
            table = table103["BASIC"];
        }
        for (size_t voltagesIndex = 0; voltagesIndex < table.size(); voltagesIndex++) {
            if (workingVoltage <= table[voltagesIndex].first) {
                clearance = table[voltagesIndex].second;
                break;
            }
        }
    }
    {
        std::vector<std::pair<double, double>> table;
        if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
            table = table104["REINFORCED"];
        }
        else {
            table = table104["BASIC"];
        }
        for (size_t voltagesIndex = 0; voltagesIndex < table.size(); voltagesIndex++) {
            if (workingVoltage <= table[voltagesIndex].first) {
                clearance = std::max(clearance, table[voltagesIndex].second);
                break;
            } 
        }
    }

    if (clearance != DBL_MAX) {
        return clearance;
    }
    else {
        throw std::invalid_argument("Too much voltage for IEC 61558-1: " + std::to_string(workingVoltage));
    }
}

double InsulationIEC61558Model::calculate_creepage_distance_over_30kHz(InsulationType insulationType, PollutionDegree pollutionDegree, double frequency, double workingVoltage) {
    std::map<double, std::vector<std::pair<double, double>>> table;

    if (workingVoltage < iec61558MinimumWorkingVoltage || insulationType == InsulationType::FUNCTIONAL) {
        return 0;
    }

    if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
        switch (pollutionDegree) {
            case PollutionDegree::P1:
                table = table108;
                break;
            case PollutionDegree::P2:
                table = table109;
                break;
            case PollutionDegree::P3:
                table = table110;
                break;
        }
    }
    else {
        switch (pollutionDegree) {
            case PollutionDegree::P1:
                table = table105;
                break;
            case PollutionDegree::P2:
                table = table106;
                break;
            case PollutionDegree::P3:
                table = table107;
                break;
        }
    }

    double previousStandardFrequency = iec61558MaximumStandardFrequency;
    std::vector<std::pair<double, double>> previousTable;
    for (auto const& [standardFrequency, voltageList] : table)
    {
        double topValue = DBL_MAX;
        if (frequency >= previousStandardFrequency && frequency <= standardFrequency) {
            topValue = linear_table_interpolation(voltageList, workingVoltage);
            if (previousTable.empty()) {
                return topValue;
            }
            else {
                auto bottomValue = linear_table_interpolation(previousTable, workingVoltage);
                double proportion = (frequency - previousStandardFrequency) / (standardFrequency - previousStandardFrequency);
                return std::lerp(bottomValue, topValue, proportion);
            }
        }
        previousStandardFrequency = standardFrequency;
        previousTable = voltageList;
    }
    throw std::invalid_argument("Too much frequency for IEC 60664-4: " + std::to_string(frequency));
}

double InsulationIEC61558Model::calculate_distance_through_insulation(InputsWrapper& inputs, bool usingThinLayers) {
    double mainSupplyVoltage = resolve_dimensional_values(inputs.get_main_supply_voltage());
    if (mainSupplyVoltage > iec61558MaximumSupplyVoltage) {
        throw std::invalid_argument("Too much supply voltage for IEC 61558-1: " + std::to_string(mainSupplyVoltage));
    }
    double workingVoltage = get_working_voltage_rms(inputs);
    auto insulationType = inputs.get_insulation_type();

    if (insulationType == InsulationType::FUNCTIONAL) {
        return 0;
    }
    double maximumFrequency = inputs.get_maximum_frequency();
    auto distanceThroughInsulation = get_distance_through_insulation_table_22(insulationType, workingVoltage, usingThinLayers);
    if (maximumFrequency > iec61558MaximumStandardFrequency) {
        distanceThroughInsulation = std::max(distanceThroughInsulation, calculate_distance_through_insulation_over_30kHz(workingVoltage));
    }

    return ceilFloat(distanceThroughInsulation, 5);
}

double InsulationIEC61558Model::calculate_withstand_voltage(InputsWrapper& inputs) {
    double mainSupplyVoltage = resolve_dimensional_values(inputs.get_main_supply_voltage());
    if (mainSupplyVoltage > iec61558MaximumSupplyVoltage) {
        throw std::invalid_argument("Too much supply voltage for IEC 61558-1: " + std::to_string(mainSupplyVoltage));
    }
    auto overvoltageCategory = inputs.get_overvoltage_category();
    double workingVoltage = get_working_voltage_rms(inputs);
    auto insulationType = inputs.get_insulation_type();
    double maximumFrequency = inputs.get_maximum_frequency();

    auto withstandVoltage = get_withstand_voltage_table_14(overvoltageCategory, insulationType, workingVoltage);

    if (maximumFrequency > iec61558MaximumStandardFrequency) {
        withstandVoltage = std::max(withstandVoltage, workingVoltage + 500);
    }

    return withstandVoltage;
}

double InsulationIEC61558Model::calculate_clearance(InputsWrapper& inputs) {
    double mainSupplyVoltage = resolve_dimensional_values(inputs.get_main_supply_voltage());
    if (mainSupplyVoltage > iec61558MaximumSupplyVoltage) {
        throw std::invalid_argument("Too much supply voltage for IEC 61558-1: " + std::to_string(mainSupplyVoltage));
    }
    auto wiringTechnology = inputs.get_design_requirements().get_wiring_technology();
    auto overvoltageCategory = inputs.get_overvoltage_category();
    double workingVoltagePeak = get_working_voltage_peak(inputs);
    double workingVoltage = get_working_voltage_rms(inputs);
    auto altitude = resolve_dimensional_values(inputs.get_altitude(), DimensionalValues::MAXIMUM);
    auto pollutionDegree = inputs.get_pollution_degree();
    auto insulationType = inputs.get_insulation_type();
    double maximumFrequency = inputs.get_maximum_frequency();
    if (insulationType == InsulationType::FUNCTIONAL) {
        return 0;
    }

    if (wiringTechnology == WiringTechnology::PRINTED) {
        pollutionDegree = PollutionDegree::P1;
    }

    auto clearance = get_clearance_table_20(overvoltageCategory, pollutionDegree, insulationType, workingVoltage);

    if (altitude > lowerAltitudeLimit) {
        if (_data.empty()) {
            auto insulationIEC60664Model = OpenMagnetics::InsulationIEC60664Model();
            return insulationIEC60664Model.calculate_clearance(inputs);
        }
        else {
            auto insulationIEC60664Model = OpenMagnetics::InsulationIEC60664Model(_data);
            return insulationIEC60664Model.calculate_clearance(inputs);
        }
    }

    if (maximumFrequency > iec61558MaximumStandardFrequency) {
        clearance = std::max(clearance, calculate_clearance_over_30kHz(insulationType, workingVoltagePeak));
    }

    return ceilFloat(clearance, 5);
}

double InsulationIEC61558Model::calculate_creepage_distance(InputsWrapper& inputs, bool includeClearance) {
    double mainSupplyVoltage = resolve_dimensional_values(inputs.get_main_supply_voltage());
    if (mainSupplyVoltage > iec61558MaximumSupplyVoltage) {
        throw std::invalid_argument("Too much supply voltage for IEC 61558-1: " + std::to_string(mainSupplyVoltage));
    }

    auto wiringTechnology = inputs.get_design_requirements().get_wiring_technology();
    double workingVoltage = get_working_voltage_rms(inputs);
    double workingVoltagePeak = get_working_voltage_peak(inputs);
    auto cti = inputs.get_cti();
    auto pollutionDegree = inputs.get_pollution_degree();
    auto insulationType = inputs.get_insulation_type();
    double maximumFrequency = inputs.get_maximum_frequency();
    if (insulationType == InsulationType::FUNCTIONAL) {
        return 0;
    }

    if (wiringTechnology == WiringTechnology::PRINTED) {
        pollutionDegree = PollutionDegree::P1;
    }

    double creepageDistance = get_creepage_distance_table_21(cti, pollutionDegree, insulationType, workingVoltage);

    if (maximumFrequency > iec61558MaximumStandardFrequency) {
        creepageDistance = std::max(creepageDistance, calculate_creepage_distance_over_30kHz(insulationType, pollutionDegree, maximumFrequency, workingVoltagePeak));
    }

    if (includeClearance) {
        creepageDistance = std::max(creepageDistance, calculate_clearance(inputs));
    }

    return ceilFloat(creepageDistance, 5);
}

double InsulationIEC60335Model::get_rated_impulse_withstand_voltage(OvervoltageCategory overvoltageCategory, double ratedVoltage){
    std::string overvoltageCategoryString = std::string{magic_enum::enum_name(overvoltageCategory)};
    auto aux = table15[overvoltageCategoryString];
    for (size_t voltagesIndex = 0; voltagesIndex < aux.size(); voltagesIndex++) {
        if (ratedVoltage <= aux[voltagesIndex].first) {
            return aux[voltagesIndex].second;
        }
    }
    throw std::invalid_argument("Too much voltage for IEC 60335-1: " + std::to_string(ratedVoltage));
}

double InsulationIEC60335Model::get_clearance_table_16(PollutionDegree pollutionDegree, WiringTechnology wiringType, InsulationType insulationType, double ratedImpulseWithstandVoltage){
    for (size_t voltagesIndex = 0; voltagesIndex < table16.size(); voltagesIndex++) {
        if (ratedImpulseWithstandVoltage <= table16[voltagesIndex].first) {
            double result = DBL_MAX;
            if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
                if (voltagesIndex < table16.size() - 1) {
                    result =  table16[voltagesIndex + 1].second;
                }
                else {
                    result =  table16[voltagesIndex].second;
                }
            }
            else {
                result =  table16[voltagesIndex].second;
            }

            if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
                if (ratedImpulseWithstandVoltage <= 800 && pollutionDegree == PollutionDegree::P3) {
                    result = 0.0008;
                }
                if (ratedImpulseWithstandVoltage <= 500 && (pollutionDegree == PollutionDegree::P1 || pollutionDegree == PollutionDegree::P2) && wiringType == WiringTechnology::PRINTED){
                    result = 0.0002;
                }
            }
            else {
                if (ratedImpulseWithstandVoltage <= 1500 && pollutionDegree == PollutionDegree::P3) {
                    result = 0.0008;
                }
                if (ratedImpulseWithstandVoltage <= 800 && (pollutionDegree == PollutionDegree::P1 || pollutionDegree == PollutionDegree::P2) && wiringType == WiringTechnology::PRINTED){
                    result = 0.0002;
                }
            }
            return result;
        }
    }
    throw std::invalid_argument("Too much voltage for IEC 60335-1: " + std::to_string(ratedImpulseWithstandVoltage));
}

double InsulationIEC60335Model::get_distance_through_insulation_table_19(OvervoltageCategory overvoltageCategory, double ratedVoltage){
    std::string overvoltageCategoryString = std::string{magic_enum::enum_name(overvoltageCategory)};
    auto aux = table19[overvoltageCategoryString];
    for (size_t voltagesIndex = 0; voltagesIndex < aux.size(); voltagesIndex++) {
        if (ratedVoltage <= aux[voltagesIndex].first) {
            return aux[voltagesIndex].second;
        }
    }
    throw std::invalid_argument("Too much voltage for IEC 60335-1: " + std::to_string(ratedVoltage));
}

double InsulationIEC60335Model::get_withstand_voltage_table_7(InsulationType insulationType, double ratedVoltage){
    std::string insulationTypeString = std::string{magic_enum::enum_name(insulationType)};
    auto aux = table7[insulationTypeString];
    for (size_t voltagesIndex = 0; voltagesIndex < aux.size(); voltagesIndex++) {
        if (ratedVoltage <= aux[voltagesIndex].first) {
            return aux[voltagesIndex].second;
        }
    }
    throw std::invalid_argument("Too much voltage for IEC 60335-1: " + std::to_string(ratedVoltage));
}

double InsulationIEC60335Model::get_withstand_voltage_formula_table_7(InsulationType insulationType, double workingVoltage){
    if (insulationType == InsulationType::BASIC) {
        return 1.2 * workingVoltage + 950;
    }
    else if (insulationType == InsulationType::SUPPLEMENTARY) {
        return 1.2 * workingVoltage + 1450;
    }
    else if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
        return 2.4 * workingVoltage + 2400;
    }
    return 0;
}


double InsulationIEC60335Model::get_creepage_distance_table_17(Cti cti, PollutionDegree pollutionDegree, double workingVoltage){
    std::string ctiString = std::string{magic_enum::enum_name(cti)};
    std::string pollutionDegreeString = std::string{magic_enum::enum_name(pollutionDegree)};
    auto aux = table17[pollutionDegreeString][ctiString];
    for (size_t voltagesIndex = 0; voltagesIndex < aux.size(); voltagesIndex++) {
        if (workingVoltage <= aux[voltagesIndex].first) {
            return aux[voltagesIndex].second;
        }
    }
    throw std::invalid_argument("Too much voltage for IEC 60335-1 Table 17: " + std::to_string(workingVoltage));
}

double InsulationIEC60335Model::get_creepage_distance_table_18(Cti cti, PollutionDegree pollutionDegree, double workingVoltage){
    std::string ctiString = std::string{magic_enum::enum_name(cti)};
    std::string pollutionDegreeString = std::string{magic_enum::enum_name(pollutionDegree)};
    auto aux = table18[pollutionDegreeString][ctiString];
    for (size_t voltagesIndex = 0; voltagesIndex < aux.size(); voltagesIndex++) {
        if (workingVoltage <= aux[voltagesIndex].first) {
            return aux[voltagesIndex].second;
        }
    }
    throw std::invalid_argument("Too much voltage for IEC 60335-1 Table 18: " + std::to_string(workingVoltage));
}


double InsulationIEC60335Model::calculate_distance_through_insulation(InputsWrapper& inputs) {
    auto insulationType = inputs.get_insulation_type();
    double mainSupplyVoltage = resolve_dimensional_values(inputs.get_main_supply_voltage());
    double maximumPrimaryVoltageRms = inputs.get_maximum_voltage_rms(0);
    auto ratedVoltage = std::max(mainSupplyVoltage, maximumPrimaryVoltageRms);
    auto overvoltageCategory = inputs.get_overvoltage_category();

    if (overvoltageCategory == OvervoltageCategory::OVC_IV) {
        throw std::invalid_argument("Overvoltage Category IV not supported in standard IEC 60335-1");
    }

    double dti = 0;
    if (insulationType == InsulationType::REINFORCED) {
        dti = std::max(0.002, get_distance_through_insulation_table_19(overvoltageCategory, ratedVoltage));
    }
    else if (insulationType == InsulationType::SUPPLEMENTARY || insulationType == InsulationType::DOUBLE) {
        dti = 0.001;
    }
 
    return dti;
}

double InsulationIEC60335Model::calculate_withstand_voltage(InputsWrapper& inputs) {
    auto insulationType = inputs.get_insulation_type();
    double mainSupplyVoltage = resolve_dimensional_values(inputs.get_main_supply_voltage());
    double maximumPrimaryVoltageRms = inputs.get_maximum_voltage_rms(0);
    double maximumVoltagePeak = inputs.get_maximum_voltage_peak();
    auto ratedVoltage = std::max(mainSupplyVoltage, maximumPrimaryVoltageRms);

    if (insulationType == InsulationType::FUNCTIONAL) {
        return 0;
    }
    else {
        double withstandVoltage = get_withstand_voltage_table_7(insulationType, ratedVoltage);
        withstandVoltage = std::max(withstandVoltage, get_withstand_voltage_formula_table_7(insulationType, maximumVoltagePeak));
        return withstandVoltage;
    }
}

double InsulationIEC60335Model::calculate_clearance(InputsWrapper& inputs) {
    auto wiringTechnology = inputs.get_design_requirements().get_wiring_technology().value();
    auto pollutionDegree = inputs.get_pollution_degree();
    double mainSupplyVoltage = resolve_dimensional_values(inputs.get_main_supply_voltage());
    double maximumPrimaryVoltageRms = inputs.get_maximum_voltage_rms(0);
    double maximumVoltageRms = inputs.get_maximum_voltage_rms();
    double maximumFrequency = inputs.get_maximum_frequency();
    auto altitude = resolve_dimensional_values(inputs.get_altitude(), DimensionalValues::MAXIMUM);
    auto ratedVoltage = std::max(mainSupplyVoltage, maximumPrimaryVoltageRms);
    auto overvoltageCategory = inputs.get_overvoltage_category();
    auto insulationType = inputs.get_insulation_type();

    if (overvoltageCategory == OvervoltageCategory::OVC_IV) {
        throw std::invalid_argument("Overvoltage Category IV not supported in standard IEC 60335-1");
    }

    double ratedImpulseWithstandVoltage = get_rated_impulse_withstand_voltage(overvoltageCategory, ratedVoltage);

    double clearance = get_clearance_table_16(pollutionDegree, wiringTechnology, insulationType, ratedImpulseWithstandVoltage);

    if (maximumVoltageRms > maximumPrimaryVoltageRms || // Boost case
        altitude > lowerAltitudeLimit ||
        maximumFrequency > ieC60335MaximumStandardFrequency) {
        double clearanceIEC60664 = DBL_MAX;
        if (_data.empty()) {
            auto insulationIEC60664Model = OpenMagnetics::InsulationIEC60664Model();
            clearanceIEC60664 = insulationIEC60664Model.calculate_clearance(inputs);
        }
        else {
            auto insulationIEC60664Model = OpenMagnetics::InsulationIEC60664Model(_data);
            clearanceIEC60664 = insulationIEC60664Model.calculate_clearance(inputs);
        }
        clearance = std::max(clearance, clearanceIEC60664);
    }

    return clearance;
}

double InsulationIEC60335Model::calculate_creepage_distance(InputsWrapper& inputs, bool includeClearance) {
    auto pollutionDegree = inputs.get_pollution_degree();
    double maximumPrimaryVoltagePeak = inputs.get_maximum_voltage_peak(0);
    double maximumFrequency = inputs.get_maximum_frequency();
    auto cti = inputs.get_cti();
    auto insulationType = inputs.get_insulation_type();

    double creepageDistance = DBL_MAX;

    if (maximumFrequency <= ieC60335MaximumStandardFrequency) {
        if (insulationType == InsulationType::FUNCTIONAL) {
            creepageDistance = get_creepage_distance_table_18(cti, pollutionDegree, maximumPrimaryVoltagePeak);
        }
        else {
            creepageDistance = get_creepage_distance_table_17(cti, pollutionDegree, maximumPrimaryVoltagePeak);
            if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
                creepageDistance *= 2;
            }
        }
    }
    else {
        double creepageDistanceIEC60664 = DBL_MAX;
        if (_data.empty()) {
            auto insulationIEC60664Model = OpenMagnetics::InsulationIEC60664Model();
            creepageDistanceIEC60664 = insulationIEC60664Model.calculate_creepage_distance(inputs);
        }
        else {
            auto insulationIEC60664Model = OpenMagnetics::InsulationIEC60664Model(_data);
            creepageDistanceIEC60664 = insulationIEC60664Model.calculate_creepage_distance(inputs);
        }
        creepageDistance = creepageDistanceIEC60664;
    }

    if (includeClearance) {
        creepageDistance = std::max(creepageDistance, calculate_clearance(inputs));
    }

    return roundFloat(creepageDistance, 5);
}

} // namespace OpenMagnetics
