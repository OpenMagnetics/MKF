#pragma once
#include <MAS.hpp>
#include "InputsWrapper.h"
#include <filesystem>
#include <fstream>
#include <iostream>


namespace OpenMagnetics {


class InsulationStandard {
  private:
  protected:
  public:
    virtual double calculate_solid_insulation(InputsWrapper inputs) = 0;
    virtual double calculate_clearance(InputsWrapper inputs) = 0;
    virtual double calculate_creepage_distance(InputsWrapper inputs, bool includeClearance = false) = 0;

    InsulationStandard() = default;
    virtual ~InsulationStandard() = default;


    // static std::shared_ptr<InsulationStandard> factory(Standard standard);
};

class InsulationCoordinator {
  private:
  protected:
  public:
    double calculate_solid_insulation(InputsWrapper inputs);
    double calculate_clearance(InputsWrapper inputs);
    double calculate_creepage_distance(InputsWrapper inputs, bool includeClearance = false);

    InsulationCoordinator() = default;
    virtual ~InsulationCoordinator() = default;

};

class InsulationIEC60664Model : public InsulationStandard {
  public:
    double calculate_solid_insulation(InputsWrapper inputs);
    double calculate_clearance(InputsWrapper inputs);
    double calculate_creepage_distance(InputsWrapper inputs, bool includeClearance = false);
    double get_rated_impulse_withstand_voltage(OvervoltageCategory overvoltageCategory, double ratedVoltage);
    double get_rated_insulation_voltage(double mainSupplyVoltage);
    double get_creepage_distance(PollutionDegree pollutionDegree, Cti cti, double voltageRms, WiringTechnology wiringType = WiringTechnology::WOUND);
    double get_creepage_distance_over_30kHz(double voltageRms, double frequency);
    std::optional<double> get_creepage_distance_planar(PollutionDegree pollutionDegree, Cti cti, double voltageRms);

    double get_clearance_table_f2(PollutionDegree pollutionDegree, InsulationType insulationType, double ratedImpulseWithstandVoltage);
    double get_clearance_table_f8(double ratedImpulseWithstandVoltage);
    std::optional<double> get_clearance_planar(double altitude, double ratedImpulseWithstandVoltage);
    double get_clearance_altitude_factor_correction(double frequency);
    double get_clearance_over_30kHz(double voltageRms, double frequency, double currentClearance);

    double iec60664Part1MaximumFrequency = 30000;
    std::vector<std::pair<double, double>> part1TableA2;
    std::map<std::string, std::vector<std::pair<double, double>>> part1TableF1;
    std::map<std::string, std::map<std::string, std::vector<std::pair<double, double>>>> part1TableF2;
    std::vector<std::pair<double, double>> part1TableF3;
    std::map<std::string, std::map<std::string, std::map<std::string, std::vector<std::pair<double, double>>>>> part1TableF5;
    std::map<std::string, std::vector<std::pair<double, double>>> part1TableF8;
    std::vector<std::pair<double, double>> part4Table1;
    std::map<double, std::vector<std::pair<double, double>>> part4Table2;
    std::map<std::string, std::vector<std::pair<double, double>>> part5Table2;
    std::map<std::string, std::vector<std::pair<double, double>>> part5Table3;
    std::map<std::string, std::map<std::string, std::vector<std::pair<double, double>>>> part5Table4;

    InsulationIEC60664Model() {
        std::string filePath = __FILE__;
        auto masPath = filePath.substr(0, filePath.rfind("/"));
        {
            auto dataFilePath = masPath + "/data/insulation_standards/IEC_60664-1.json";
            std::ifstream jsonFile(dataFilePath);
            json jf = json::parse(jsonFile);
            part1TableA2 = jf["A.2"];
            part1TableF1 = jf["F.1"];
            part1TableF2 = jf["F.2"];
            part1TableF3 = jf["F.3"];
            part1TableF5 = jf["F.5"];
            part1TableF8 = jf["F.8"];
        }
        {
            auto dataFilePath = masPath + "/data/insulation_standards/IEC_60664-4.json";
            std::ifstream jsonFile(dataFilePath);
            json jf = json::parse(jsonFile);
            part4Table1 = jf["Table 1"];

            std::map<std::string, std::vector<std::pair<double, double>>> temp;
            temp = jf["Table 2"];
            for (auto const& [standardFrequencyString, voltageList] : temp)
            {
                double standardFrequency = std::stod(standardFrequencyString);
                part4Table2[standardFrequency] = voltageList;
            }
        }
        {
            auto dataFilePath = masPath + "/data/insulation_standards/IEC_60664-5.json";
            std::ifstream jsonFile(dataFilePath);
            json jf = json::parse(jsonFile);
            part5Table2 = jf["Table 2"];
            part5Table3 = jf["Table 3"];
            part5Table4 = jf["Table 4"];
        }
    }

    InsulationIEC60664Model(json data) {
        part1TableA2 = data["IEC_60664-1"]["A.2"];
        part1TableF1 = data["IEC_60664-1"]["F.1"];
        part1TableF2 = data["IEC_60664-1"]["F.2"];
        part1TableF3 = data["IEC_60664-1"]["F.3"];
        part1TableF5 = data["IEC_60664-1"]["F.5"];
        part1TableF8 = data["IEC_60664-1"]["F.8"];

        part4Table1 = data["IEC_60664-4"]["Table 1"];
        std::map<std::string, std::vector<std::pair<double, double>>> temp;
        temp = data["IEC_60664-4"]["Table 2"];
        for (auto const& [standardFrequencyString, voltageList] : temp)
        {
            double standardFrequency = std::stod(standardFrequencyString);
            part4Table2[standardFrequency] = voltageList;
        }

        part5Table2 = data["IEC_60664-5"]["Table 2"];
        part5Table3 = data["IEC_60664-5"]["Table 3"];
        part5Table4 = data["IEC_60664-5"]["Table 4"];
    }
};

class InsulationIEC62368Model : public InsulationStandard {
  public:
    double calculate_solid_insulation(InputsWrapper inputs);
    double calculate_clearance(InputsWrapper inputs);
    double calculate_creepage_distance(InputsWrapper inputs, bool includeClearance = false);

    double get_working_voltage(InputsWrapper inputs);
    double get_working_voltage_rms(InputsWrapper inputs);
    double get_required_withstand_voltage(InputsWrapper inputs);
    double get_voltage_due_to_transient_overvoltages(double requiredWithstandVoltage, InsulationType insulationType);
    double get_voltage_due_to_recurring_peak_voltages(double workingVoltage, InsulationType insulationType);
    double get_voltage_due_to_temporary_overvoltages(double supplyVoltageRms, InsulationType insulationType);
    double get_reduction_factor_per_material(std::string material, double frequency);
    double get_clearance_table_10(double supplyVoltagePeak, InsulationType insulationType, PollutionDegree pollutionDegree);
    double get_clearance_table_11(double supplyVoltagePeak, InsulationType insulationType, PollutionDegree pollutionDegree);
    double get_clearance_table_14(double supplyVoltagePeak, InsulationType insulationType, PollutionDegree pollutionDegree);
    double get_altitude_factor(double altitude);
    double get_creepage_distance_table_17(double voltageRms, InsulationType insulationType, PollutionDegree pollutionDegree, Cti cti);
    double get_creepage_distance_table_18(double voltageRms, double frequency, PollutionDegree pollutionDegree, InsulationType insulationType);
    double get_voltage_due_to_temporary_overvoltages_procedure_1(double supplyVoltagePeak);
    double get_mains_transient_voltage(double supplyVoltagePeak, OvervoltageCategory overvoltageCategory);
    double get_distance_table_G13(double workingVoltage, InsulationType insulationType);

    std::map<std::string, std::map<std::string, std::vector<std::pair<double, double>>>> table10;
    std::map<std::string, std::vector<std::pair<double, double>>> table11;
    std::map<std::string, std::vector<std::pair<double, double>>> table12;
    std::map<std::string, std::map<std::string, std::vector<std::pair<double, double>>>> table14;
    std::vector<std::pair<double, double>> table15;
    std::vector<std::pair<double, double>> table16;
    std::map<std::string, std::map<std::string, std::vector<std::pair<double, double>>>> table17;
    std::map<std::string, std::vector<std::pair<double, double>>> table18;
    std::map<std::string, std::vector<std::pair<double, double>>> table21;
    std::map<std::string, std::vector<std::pair<double, double>>> table22;
    std::map<std::string, std::vector<std::pair<double, double>>> table25;
    std::map<std::string, std::vector<std::pair<double, double>>> table26;
    std::map<std::string, std::vector<std::pair<double, double>>> table27;
    std::map<std::string, std::vector<std::pair<double, double>>> tableG4;
    std::map<std::string, std::vector<std::pair<double, double>>> tableG13;

    InsulationIEC62368Model() {
        std::string filePath = __FILE__;
        auto masPath = filePath.substr(0, filePath.rfind("/"));
        {
            auto dataFilePath = masPath + "/data/insulation_standards/IEC_62368-1.json";
            std::ifstream jsonFile(dataFilePath);
            json jf = json::parse(jsonFile);
            table10 = jf["Table 10"];
            table11 = jf["Table 11"];
            table12 = jf["Table 12"];
            table14 = jf["Table 14"];
            table15 = jf["Table 15"];
            table16 = jf["Table 16"];
            table17 = jf["Table 17"];
            table18 = jf["Table 18"];
            table21 = jf["Table 21"];
            table22 = jf["Table 22"];
            table25 = jf["Table 25"];
            table26 = jf["Table 26"];
            table27 = jf["Table 27"];
            tableG4 = jf["Table G.4"];
            tableG13 = jf["Table G.13"];
        }
    }

    InsulationIEC62368Model(json data) {
        table10 = data["Table 10"];
        table11 = data["Table 11"];
        table12 = data["Table 12"];
        table14 = data["Table 14"];
        table15 = data["Table 15"];
        table16 = data["Table 16"];
        table17 = data["Table 17"];
        table18 = data["Table 18"];
        table21 = data["Table 21"];
        table25 = data["Table 25"];
        table26 = data["Table 26"];
        table27 = data["Table 27"];
        tableG4 = data["Table G.4"];
        tableG13 = data["Table G.13"];
    }
};

} // namespace OpenMagnetics