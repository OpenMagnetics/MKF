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
    virtual double calculate_creepage_distance(InputsWrapper inputs) = 0;

    InsulationStandard() = default;
    virtual ~InsulationStandard() = default;


    // static std::shared_ptr<InsulationStandard> factory(Standard standard);
};

class InsulationIEC60664Model : public InsulationStandard {
  public:
    double calculate_solid_insulation(InputsWrapper inputs);
    double calculate_clearance(InputsWrapper inputs);
    double calculate_creepage_distance(InputsWrapper inputs);
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
};

} // namespace OpenMagnetics