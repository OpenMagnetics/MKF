#include "Utils.h"
#include "Insulation.h"
#include <magic_enum.hpp>
#include "spline.h"
#include "WindingSkinEffectLosses.h"
#include <cfloat>


namespace OpenMagnetics {


double InsulationCoordinator::calculate_solid_insulation(InputsWrapper inputs) {
    double solidInsulation = 0;
    for (auto standard : inputs.get_standards()) {
        switch (standard) {
            case InsulationStandards::IEC_603351: {
                throw std::invalid_argument("IEC_60335 not implemented yet");
                break;
            }
            case InsulationStandards::IEC_606641: {
                auto standardObj = OpenMagnetics::InsulationIEC60664Model();
                solidInsulation = std::max(solidInsulation, standardObj.calculate_solid_insulation(inputs));
                break;
            }
            case InsulationStandards::IEC_615581: {
                throw std::invalid_argument("IEC_61558 not implemented yet");
                break;
            }
            case InsulationStandards::IEC_623681: {
                auto standardObj = OpenMagnetics::InsulationIEC62368Model();
                solidInsulation = std::max(solidInsulation, standardObj.calculate_solid_insulation(inputs));
                break;
            }
        }
    }
    return solidInsulation;
}

double InsulationCoordinator::calculate_clearance(InputsWrapper inputs) {
    double clearance = 0;
    for (auto standard : inputs.get_standards()) {
        switch (standard) {
            case InsulationStandards::IEC_603351: {
                throw std::invalid_argument("IEC_60335 not implemented yet");
                break;
            }
            case InsulationStandards::IEC_606641: {
                auto standardObj = OpenMagnetics::InsulationIEC60664Model();
                clearance = std::max(clearance, standardObj.calculate_clearance(inputs));
                break;
            }
            case InsulationStandards::IEC_615581: {
                throw std::invalid_argument("IEC_61558 not implemented yet");
                break;
            }
            case InsulationStandards::IEC_623681: {
                auto standardObj = OpenMagnetics::InsulationIEC62368Model();
                clearance = std::max(clearance, standardObj.calculate_clearance(inputs));
                break;
            }
        }
    }
    return clearance;
}

double InsulationCoordinator::calculate_creepage_distance(InputsWrapper inputs, bool includeClearance) {
    double creepageDistance = 0;
    for (auto standard : inputs.get_standards()) {
        switch (standard) {
            case InsulationStandards::IEC_603351: {
                throw std::invalid_argument("IEC_60335 not implemented yet");
                break;
            }
            case InsulationStandards::IEC_606641: {
                auto standard = OpenMagnetics::InsulationIEC60664Model();
                creepageDistance = std::max(creepageDistance, standard.calculate_creepage_distance(inputs, includeClearance));
                break;
            }
            case InsulationStandards::IEC_615581: {
                throw std::invalid_argument("IEC_61558 not implemented yet");
                break;
            }
            case InsulationStandards::IEC_623681: {
                auto standard = OpenMagnetics::InsulationIEC62368Model();
                creepageDistance = std::max(creepageDistance, standard.calculate_creepage_distance(inputs, includeClearance));
                break;
            }
        }
    }
    return creepageDistance;
}

double InsulationIEC60664Model::get_rated_impulse_withstand_voltage(OvervoltageCategory overvoltageCategory, double ratedVoltage){
    std::string overvoltageCategoryString = std::string{magic_enum::enum_name(overvoltageCategory)};
    auto aux = part1TableF1[overvoltageCategoryString];
    for (size_t voltagesIndex = 0; voltagesIndex < aux.size(); voltagesIndex++) {
        if (ratedVoltage <= aux[voltagesIndex].first) {
            return aux[voltagesIndex].second;
        }
    }
    throw std::invalid_argument("Too much voltage for IEC 60664-1: " + std::to_string(ratedVoltage));
}

double InsulationIEC60664Model::get_clearance_table_f2(PollutionDegree pollutionDegree, InsulationType insulationType, double ratedImpulseWithstandVoltage){
    std::string pollutionDegreeString = std::string{magic_enum::enum_name(pollutionDegree)};
    auto aux = part1TableF2["inhomogeneusField"][pollutionDegreeString];
    for (size_t voltagesIndex = 0; voltagesIndex < aux.size(); voltagesIndex++) {
        if (ratedImpulseWithstandVoltage <= aux[voltagesIndex].first) {
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
    std::vector<std::pair<double, double>> aux;
    if (altitude <= 2000) {
        aux = part5Table2["inhomogeneusField"];
    }
    else {
        aux = part5Table3["inhomogeneusField"];
    }

    std::vector<double> x, y;
    bool insideTable = false;
    for (auto& voltagePair : aux) {
        if (ratedImpulseWithstandVoltage < voltagePair.first) {
            insideTable = true;
        }
        x.push_back(voltagePair.first);
        y.push_back(voltagePair.second);
    }

    if (insideTable) {
        tk::spline interp(x, y);
        return interp(ratedImpulseWithstandVoltage);
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
    for (auto const& [standardFrequency, voltageList] : part4Table2)
    {
        if (frequency >= previousStandardFrequency && frequency <= standardFrequency) {
            int n = voltageList.size();
            std::vector<double> x, y;

            for (int i = 0; i < n; i++) {
                x.push_back(voltageList[i].first);
                y.push_back(voltageList[i].second);
            }

            tk::spline interp(x, y);
            return interp(voltageRms);
        }
        previousStandardFrequency = standardFrequency;
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

    auto aux = part5Table4[pollutionDegreeString][ctiString];

    std::vector<double> x, y;
    bool insideTable = false;
    for (auto& voltagePair : aux) {
        if (voltageRms < voltagePair.first) {
            insideTable = true;
        }
        x.push_back(voltagePair.first);
        y.push_back(voltagePair.second);
    }

    if (insideTable) {
        tk::spline interp(x, y);
        return interp(voltageRms);
    }

    // If we find to option in this table, we go back to IEC 60664-1
    return std::nullopt;
}

double InsulationIEC60664Model::get_clearance_altitude_factor_correction(double altitude){
    int n = part1TableA2.size();
    std::vector<double> x, y;

    for (int i = 0; i < n; i++) {
        x.push_back(part1TableA2[i].first);
        y.push_back(part1TableA2[i].second);
    }

    tk::spline interp(x, y);
    return interp(altitude);
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
        int n = part4Table1.size();
        std::vector<double> x, y;

        for (int i = 0; i < n; i++) {
            x.push_back(part4Table1[i].first);
            y.push_back(part4Table1[i].second);
        }

        tk::spline interp(x, y);
        return interp(ratedVoltagePeak);
    }
    return currentClearance;
}

double InsulationIEC60664Model::calculate_solid_insulation(InputsWrapper inputs) {
    double maximumVoltageRms = inputs.get_maximum_voltage_rms();
    auto overvoltageCategory = inputs.get_overvoltage_category();
    auto insulationType = inputs.get_insulation_type();

    double voltageDueToTransientOvervoltages = get_rated_impulse_withstand_voltage(overvoltageCategory, maximumVoltageRms);
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

    // TODO take into account greater than 30kHzz case
    return std::max(voltageDueToTransientOvervoltages, std::max(voltageDueToTemporaryWithstandOvervoltages, std::max(voltageDueToRecurringPeakVoltages, voltageDueToSteadyStateVoltages)));
}

double InsulationIEC60664Model::calculate_clearance(InputsWrapper inputs) {
    auto wiringTechnology = inputs.get_design_requirements().get_wiring_technology();
    auto pollutionDegree = inputs.get_pollution_degree();
    double maximumFrequency = inputs.get_maximum_frequency();
    double ratedVoltage = inputs.get_maximum_voltage_rms();
    auto overvoltageCategory = inputs.get_overvoltage_category();
    auto altitude = resolve_dimensional_values(inputs.get_altitude(), DimensionalValues::MAXIMUM);
    auto insulationType = inputs.get_insulation_type();
    double ratedImpulseWithstandVoltage = get_rated_impulse_withstand_voltage(overvoltageCategory, ratedVoltage);
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
        clearanceToWithstandTransientOvervoltages = get_clearance_table_f2(pollutionDegree, insulationType, ratedImpulseWithstandVoltage);
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

    if (altitude > 2000) {
        clearance *= get_clearance_altitude_factor_correction(altitude);
    }

    return clearance;
}

double InsulationIEC60664Model::calculate_creepage_distance(InputsWrapper inputs, bool includeClearance) {
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
            double creepageDistanceOver30kHz = get_creepage_distance_over_30kHz(voltageRms, maximumFrequency);
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
        auto ratedImpulseWithstandVoltage = get_rated_impulse_withstand_voltage(overvoltageCategory, ratedVoltage);
        auto clearanceToWithstandTransientOvervoltages = get_clearance_table_f2(pollutionDegree, insulationType, ratedImpulseWithstandVoltage);
        if (clearanceToWithstandTransientOvervoltages < creepageDistance) {
            allowSmallerCreepageDistanceThanClearance = true;
        }
    }

    if (!allowSmallerCreepageDistanceThanClearance && includeClearance) {
        creepageDistance = std::max(creepageDistance, calculate_clearance(inputs));
    }

    return roundFloat(creepageDistance, 5);
}

double InsulationIEC62368Model::get_working_voltage(InputsWrapper inputs) {
    return inputs.get_maximum_voltage_peak();
}

double InsulationIEC62368Model::get_working_voltage_rms(InputsWrapper inputs) {
    return inputs.get_maximum_voltage_rms();
}

double InsulationIEC62368Model::get_required_withstand_voltage(InputsWrapper inputs) {
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
        table = table25["Reinforced"];
    }
    else {
        table = table25["Basic"];
    }
    for (auto& voltagePair : table) {
        if (requiredWithstandVoltage < voltagePair.first) {
            return voltagePair.second;
        }
    }
    throw std::invalid_argument("Too much voltage for IEC 62368-1: " + std::to_string(requiredWithstandVoltage));
}

double InsulationIEC62368Model::get_voltage_due_to_recurring_peak_voltages(double workingVoltage, InsulationType insulationType){
    std::vector<std::pair<double, double>> table;
    if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
        table = table26["Reinforced"];
    }
    else {
        table = table26["Basic"];
    }
    for (auto& voltagePair : table) {
        if (workingVoltage < voltagePair.first) {
            return voltagePair.second;
        }
    }
    throw std::invalid_argument("Too much voltage for IEC 62368-1: " + std::to_string(workingVoltage));
}

double InsulationIEC62368Model::get_voltage_due_to_temporary_overvoltages(double supplyVoltageRms, InsulationType insulationType){
    std::vector<std::pair<double, double>> table;
    if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
        table = table27["Reinforced"];
    }
    else {
        table = table27["Basic"];
    }
    for (auto& voltagePair : table) {
        if (supplyVoltageRms < voltagePair.first) {
            return voltagePair.second;
        }
    }
    throw std::invalid_argument("Too much voltage for IEC 62368-1: " + std::to_string(supplyVoltageRms));
}

double InsulationIEC62368Model::get_reduction_factor_per_material(std::string material, double frequency){
    double previousStandardFrequency = 30000;
    for (auto const& [standardFrequency, reductionFactor] : table22[material])
    {
        if (frequency >= previousStandardFrequency && frequency <= standardFrequency) {
            return reductionFactor;
        }
        previousStandardFrequency = standardFrequency;
    }
    throw std::invalid_argument("Too much frequency for IEC 62368-1: " + std::to_string(frequency));
}

double InsulationIEC62368Model::get_clearance_table_10(double supplyVoltagePeak, InsulationType insulationType, PollutionDegree pollutionDegree){
    std::string pollutionDegreeString = std::string{magic_enum::enum_name(pollutionDegree)};
    std::vector<std::pair<double, double>> table;
    if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
        table = table10["Reinforced"][pollutionDegreeString];
    }
    else {
        table = table10["Basic"][pollutionDegreeString];
    }

    std::vector<double> x, y;
    for (auto& voltagePair : table) {
        x.push_back(voltagePair.first);
        y.push_back(voltagePair.second);
    }

    tk::spline interp(x, y);
    double result = interp(supplyVoltagePeak);

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
        table = table14["Reinforced"][pollutionDegreeString];
    }
    else {
        table = table14["Basic"][pollutionDegreeString];
    }

    std::vector<double> x, y;
    for (auto& voltagePair : table) {
        x.push_back(voltagePair.first);
        y.push_back(voltagePair.second);
    }

    tk::spline interp(x, y);
    double result = interp(supplyVoltagePeak);
    
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
        table = table11["Reinforced"];
    }
    else {
        table = table11["Basic"];
    }

    std::vector<double> x, y;
    for (auto& voltagePair : table) {
        x.push_back(voltagePair.first);
        y.push_back(voltagePair.second);
    }

    tk::spline interp(x, y);
    double valuePollutionDegree2 = interp(supplyVoltagePeak);
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
        table = tableG13["Reinforced"];
    }
    else {
        table = tableG13["Basic"];
    }

    std::vector<double> x, y;
    for (auto& voltagePair : table) {
        x.push_back(voltagePair.first);
        y.push_back(voltagePair.second);
    }

    tk::spline interp(x, y);
    double result = interp(workingVoltage);
    return ceilFloat(result, 4);
}


double InsulationIEC62368Model::get_creepage_distance_table_17(double voltageRms, InsulationType insulationType, PollutionDegree pollutionDegree, Cti cti){
    std::string pollutionDegreeString = std::string{magic_enum::enum_name(pollutionDegree)};
    std::string ctiString = std::string{magic_enum::enum_name(cti)};
    std::vector<std::pair<double, double>> table = table17[pollutionDegreeString][ctiString];

    std::vector<double> x, y;
    for (auto& voltagePair : table) {
        x.push_back(voltagePair.first);
        y.push_back(voltagePair.second);
    }

    tk::spline interp(x, y);
    double valueBasic = interp(voltageRms);

    if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
        return ceilFloat(valueBasic * 2, 4);
    }
    else {
        return ceilFloat(valueBasic, 4);
    }
}


double InsulationIEC62368Model::get_creepage_distance_table_18(double voltageRms, double frequency, PollutionDegree pollutionDegree, InsulationType insulationType){
    double previousStandardFrequency = 30000;
    double valuePollutionDegree1 = DBL_MAX;
    for (auto const& [standardFrequencyString, voltageList] : table18)
    {
        double standardFrequency = std::stod(standardFrequencyString);
        if ((frequency >= previousStandardFrequency && frequency <= standardFrequency) || standardFrequency == 400000) {
            int n = voltageList.size();
            std::vector<double> x, y;

            for (int i = 0; i < n; i++) {
                x.push_back(voltageList[i].first);
                y.push_back(voltageList[i].second);
            }

            tk::spline interp(x, y);
            valuePollutionDegree1 = interp(voltageRms);
            break;
        }
        previousStandardFrequency = standardFrequency;
    }
    if (valuePollutionDegree1 == DBL_MAX) {
        throw std::invalid_argument("Too much frequency for IEC 62368-1: " + std::to_string(frequency));
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
    std::vector<double> x, y;
    for (auto& voltagePair : table16) {
        x.push_back(voltagePair.first);
        y.push_back(voltagePair.second);
    }

    tk::spline interp(x, y);
    return ceilFloat(interp(altitude), 5);
}

double InsulationIEC62368Model::get_mains_transient_voltage(double supplyVoltagePeak, OvervoltageCategory overvoltageCategory){
    std::string overvoltageCategoryString = std::string{magic_enum::enum_name(overvoltageCategory)};
    auto table = table12[overvoltageCategoryString];
    for (size_t voltagesIndex = 0; voltagesIndex < table.size(); voltagesIndex++) {
        if (supplyVoltagePeak <= table[voltagesIndex].first) {
            return table[voltagesIndex].second;
        }
    }
    throw std::invalid_argument("Too much voltage for IEC 62368-1: " + std::to_string(supplyVoltagePeak));
}

double InsulationIEC62368Model::calculate_solid_insulation(InputsWrapper inputs) {
    double maximumFrequency = inputs.get_maximum_frequency();
    double workingVoltage = get_working_voltage(inputs);
    double requiredWithstandVoltage = get_required_withstand_voltage(inputs);
    double mainSupplyVoltage = resolve_dimensional_values(inputs.get_main_supply_voltage());
    auto insulationType = inputs.get_insulation_type();

    if (maximumFrequency > 30000) {
        auto reductionFactor = get_reduction_factor_per_material("Other thin foil materials", maximumFrequency);
        if (insulationType == InsulationType::REINFORCED || insulationType == InsulationType::DOUBLE) {
            mainSupplyVoltage = 1.2 * 2 * mainSupplyVoltage / reductionFactor;
        }
        else {
            mainSupplyVoltage = 1.2 * mainSupplyVoltage / reductionFactor;
        }
    }

    double voltageDueToTransientOvervoltages = get_voltage_due_to_transient_overvoltages(workingVoltage, insulationType);
    double voltageDueToRecurringPeakVoltages = get_voltage_due_to_recurring_peak_voltages(requiredWithstandVoltage, insulationType);
    double voltageDueToTemporaryOvervoltages = get_voltage_due_to_temporary_overvoltages(mainSupplyVoltage, insulationType);

    return std::max(voltageDueToTransientOvervoltages, std::max(voltageDueToRecurringPeakVoltages, voltageDueToTemporaryOvervoltages));
}

double InsulationIEC62368Model::calculate_clearance(InputsWrapper inputs) {
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
    if (maximumFrequency <= 30000) {
        clearanceProcedure1 = get_clearance_table_10(voltageProcedure1, insulationType, pollutionDegree);
    }
    else {
        clearanceProcedure1 = get_clearance_table_11(voltagePeak, insulationType, pollutionDegree);
    }
    double mainsTransientVoltage = get_mains_transient_voltage(mainSupplyVoltage, overvoltageCategory);
    clearanceProcedure2 = get_clearance_table_14(std::max(requiredWithstandVoltage, mainsTransientVoltage), insulationType, pollutionDegree);

    double altitudeFactor = get_altitude_factor(altitude);

    return ceilFloat(altitudeFactor * std::max(clearanceProcedure1, clearanceProcedure2), 5);
    // TODO: remove distance if FIW complies woith conditions in table G.4
}

double InsulationIEC62368Model::calculate_creepage_distance(InputsWrapper inputs, bool includeClearance) {
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

    if (maximumFrequency <= 30000) {
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

    // TODO: remove distance if FIW complies woith conditions in table G.4
}

} // namespace OpenMagnetics
