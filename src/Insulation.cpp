#include "Utils.h"
#include "Insulation.h"
#include <magic_enum.hpp>
#include "spline.h"
#include "WindingSkinEffectLosses.h"

namespace OpenMagnetics {


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
        tk::spline interp(x, y, tk::spline::cspline_hermite);
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
        if (frequency > previousStandardFrequency && frequency <= standardFrequency) {
            int n = voltageList.size();
            std::vector<double> x, y;

            for (int i = 0; i < n; i++) {
                x.push_back(voltageList[i].first);
                y.push_back(voltageList[i].second);
            }

            tk::spline interp(x, y, tk::spline::cspline_hermite);
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
        tk::spline interp(x, y, tk::spline::cspline_hermite);
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

    tk::spline interp(x, y, tk::spline::cspline_hermite);
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

        tk::spline interp(x, y, tk::spline::cspline_hermite);
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

    // std::cout << "ratedVoltage:" << ratedVoltage << std::endl;
    // std::cout << "ratedImpulseWithstandVoltage:" << ratedImpulseWithstandVoltage << std::endl;
    // std::cout << "clearanceToWithstandTransientOvervoltages:" << clearanceToWithstandTransientOvervoltages << std::endl;
    // std::cout << "clearanceToWithstandSteadyStatePeakVoltages:" << clearanceToWithstandSteadyStatePeakVoltages << std::endl;

    double clearance = std::max(clearanceToWithstandTransientOvervoltages, clearanceToWithstandSteadyStatePeakVoltages);

    if (altitude > 2000) {
        clearance *= get_clearance_altitude_factor_correction(altitude);
    }

    return clearance;
}

double InsulationIEC60664Model::calculate_creepage_distance(InputsWrapper inputs) {
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
    return roundFloat(creepageDistance, 5);
}

} // namespace OpenMagnetics
