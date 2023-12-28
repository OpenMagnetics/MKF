#pragma once
#include <MAS.hpp>
#include "InputsWrapper.h"
#include "WireWrapper.h"
#include "InsulationMaterialWrapper.h"
#include <filesystem>
#include <fstream>
#include <iostream>


namespace OpenMagnetics {


class CoilSectionInterface {
    public:
    enum class LayerPurpose : int { MECHANICAL, INSULATING };
    CoilSectionInterface() = default;
    virtual ~CoilSectionInterface() = default;

    private:
    double totalMarginTapeDistance;
    double solidInsulationThickness;
    size_t numberLayersInsulation;
    LayerPurpose layerPurpose;

    public:

    double get_total_margin_tape_distance() const { return totalMarginTapeDistance; }
    void set_total_margin_tape_distance(const double & value) { this->totalMarginTapeDistance = value; }

    double get_solid_insulation_thickness() const { return solidInsulationThickness; }
    void set_solid_insulation_thickness(const double & value) { this->solidInsulationThickness = value; }

    size_t get_number_layers_insulation() const { return numberLayersInsulation; }
    void set_number_layers_insulation(const size_t & value) { this->numberLayersInsulation = value; }

    LayerPurpose get_layer_purpose() const { return layerPurpose; }
    void set_layer_purpose(LayerPurpose value) { this->layerPurpose = value; }
};

class InsulationStandard {
  private:
  protected:
  public:
    virtual double calculate_withstand_voltage(InputsWrapper inputs) = 0;
    virtual double calculate_clearance(InputsWrapper inputs) = 0;
    virtual double calculate_creepage_distance(InputsWrapper inputs, bool includeClearance = false) = 0;

    InsulationStandard() = default;
    virtual ~InsulationStandard() = default;
    double iec60664Part1MaximumFrequency = 30000;
    double lowerAltitudeLimit = 2000;


    // static std::shared_ptr<InsulationStandard> factory(Standard standard);
};

class InsulationIEC60664Model : public InsulationStandard {
  public:
    double calculate_distance_through_insulation(InputsWrapper inputs);
    double calculate_withstand_voltage(InputsWrapper inputs);
    double calculate_clearance(InputsWrapper inputs);
    double calculate_creepage_distance(InputsWrapper inputs, bool includeClearance = false);
    double get_rated_impulse_withstand_voltage(OvervoltageCategory overvoltageCategory, double ratedVoltage, InsulationType insulationType);
    double get_rated_insulation_voltage(double mainSupplyVoltage);
    double get_creepage_distance(PollutionDegree pollutionDegree, Cti cti, double voltageRms, WiringTechnology wiringType = WiringTechnology::WOUND);
    double get_creepage_distance_over_30kHz(double voltageRms, double frequency);
    std::optional<double> get_creepage_distance_planar(PollutionDegree pollutionDegree, Cti cti, double voltageRms);

    double get_clearance_table_f2(PollutionDegree pollutionDegree, double ratedImpulseWithstandVoltage);
    double get_clearance_table_f8(double ratedImpulseWithstandVoltage);
    std::optional<double> get_clearance_planar(double altitude, double ratedImpulseWithstandVoltage);
    double get_clearance_altitude_factor_correction(double frequency);
    double get_clearance_over_30kHz(double voltageRms, double frequency, double currentClearance);
    double calculate_distance_through_insulation_over_30kHz(double workingVoltage);
    bool electric_field_strenth_is_valid(double dti, double voltage);

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
    double calculate_withstand_voltage(InputsWrapper inputs);
    double calculate_clearance(InputsWrapper inputs);
    double calculate_creepage_distance(InputsWrapper inputs, bool includeClearance = false);
    double calculate_distance_through_insulation(InputsWrapper inputs);

    double get_es2_voltage_limit(double frequency);
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

    double iec62368LowerFrequency = 30000;
    std::map<std::string, std::map<std::string, std::vector<std::pair<double, double>>>> table10;
    std::map<std::string, std::vector<std::pair<double, double>>> table11;
    std::map<std::string, std::vector<std::pair<double, double>>> table12;
    std::map<std::string, std::map<std::string, std::vector<std::pair<double, double>>>> table14;
    std::vector<std::pair<double, double>> table15;
    std::vector<std::pair<double, double>> table16;
    std::map<std::string, std::map<std::string, std::vector<std::pair<double, double>>>> table17;
    std::map<double, std::vector<std::pair<double, double>>> table18;
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

            std::map<std::string, std::vector<std::pair<double, double>>> temp;
            temp = jf["Table 18"];
            for (auto const& [standardFrequencyString, voltageList] : temp)
            {
                double standardFrequency = std::stod(standardFrequencyString);
                table18[standardFrequency] = voltageList;
            }

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
        table10 = data["IEC_62368-1"]["Table 10"];
        table11 = data["IEC_62368-1"]["Table 11"];
        table12 = data["IEC_62368-1"]["Table 12"];
        table14 = data["IEC_62368-1"]["Table 14"];
        table15 = data["IEC_62368-1"]["Table 15"];
        table16 = data["IEC_62368-1"]["Table 16"];
        table17 = data["IEC_62368-1"]["Table 17"];

        std::map<std::string, std::vector<std::pair<double, double>>> temp;
        temp = data["IEC_62368-1"]["Table 18"];
        for (auto const& [standardFrequencyString, voltageList] : temp)
        {
            double standardFrequency = std::stod(standardFrequencyString);
            table18[standardFrequency] = voltageList;
        }

        table21 = data["IEC_62368-1"]["Table 21"];
        table22 = data["IEC_62368-1"]["Table 22"];
        table25 = data["IEC_62368-1"]["Table 25"];
        table26 = data["IEC_62368-1"]["Table 26"];
        table27 = data["IEC_62368-1"]["Table 27"];
        tableG4 = data["IEC_62368-1"]["Table G.4"];
        tableG13 = data["IEC_62368-1"]["Table G.13"];
    }
};

class InsulationIEC61558Model : public InsulationStandard {
  private:
    json _data;
  public:
    double calculate_distance_through_insulation(InputsWrapper inputs, bool usingThinLayers = true);
    double calculate_withstand_voltage(InputsWrapper inputs);
    double calculate_clearance(InputsWrapper inputs);
    double calculate_creepage_distance(InputsWrapper inputs, bool includeClearance = false);
    double get_withstand_voltage_table_14(OvervoltageCategory overvoltageCategory, InsulationType insulationType, double workingVoltage);
    double get_clearance_table_20(OvervoltageCategory overvoltageCategory, PollutionDegree pollutionDegree, InsulationType insulationType, double workingVoltage);
    double get_creepage_distance_table_21(Cti cti, PollutionDegree pollutionDegree, InsulationType insulationType, double workingVoltage);
    double get_distance_through_insulation_table_22(InsulationType insulationType, double workingVoltage, bool usingThinLayers = true);
    double get_working_voltage_rms(InputsWrapper inputs);
    double get_working_voltage_peak(InputsWrapper inputs);
    double calculate_distance_through_insulation_over_30kHz(double workingVoltage);
    double calculate_clearance_over_30kHz(InsulationType insulationType, double workingVoltage);
    double calculate_creepage_distance_over_30kHz(InsulationType insulationType, PollutionDegree pollutionDegree, double frequency, double workingVoltage);
    bool electric_field_strenth_is_valid(double dti, double voltage);

    double iec61558MinimumWorkingVoltage = 25;
    double iec61558MaximumSupplyVoltage = 1100;
    double iec61558MaximumStandardFrequency = 30000;
    std::map<std::string, std::map<std::string, std::vector<std::pair<double, double>>>> table14;
    std::map<std::string, std::map<std::string, std::map<std::string, std::vector<std::pair<double, double>>>>> table20;
    std::map<std::string, std::map<std::string, std::map<std::string, std::vector<std::pair<double, double>>>>> table21;
    std::map<std::string, std::map<std::string, std::vector<std::pair<double, double>>>> table22;
    std::map<std::string, std::vector<std::pair<double, double>>> table23;
    std::map<std::string, std::vector<std::pair<double, double>>> table102;
    std::map<std::string, std::vector<std::pair<double, double>>> table103;
    std::map<std::string, std::vector<std::pair<double, double>>> table104;
    std::map<double, std::vector<std::pair<double, double>>> table105;
    std::map<double, std::vector<std::pair<double, double>>> table106;
    std::map<double, std::vector<std::pair<double, double>>> table107;
    std::map<double, std::vector<std::pair<double, double>>> table108;
    std::map<double, std::vector<std::pair<double, double>>> table109;
    std::map<double, std::vector<std::pair<double, double>>> table110;


    InsulationIEC61558Model() {
        std::string filePath = __FILE__;
        auto masPath = filePath.substr(0, filePath.rfind("/"));
        {
            auto dataFilePath = masPath + "/data/insulation_standards/IEC_61558-1.json";
            std::ifstream jsonFile(dataFilePath);
            json jf = json::parse(jsonFile);
            table14 = jf["Table 14"];
            table20 = jf["Table 20"];
            table21 = jf["Table 21"];
            table22 = jf["Table 22"];
            table23 = jf["Table 23"];
        }
        {
            auto dataFilePath = masPath + "/data/insulation_standards/IEC_61558-2-16.json";
            std::ifstream jsonFile(dataFilePath);
            json jf = json::parse(jsonFile);
            table102 = jf["Table 102"];
            table103 = jf["Table 103"];
            table104 = jf["Table 104"];


            std::map<std::string, std::vector<std::pair<double, double>>> temp;
            temp = jf["Table 105"];
            for (auto const& [standardFrequencyString, voltageList] : temp)
            {
                double standardFrequency = std::stod(standardFrequencyString);
                table105[standardFrequency] = voltageList;
            }
            temp = jf["Table 106"];
            for (auto const& [standardFrequencyString, voltageList] : temp)
            {
                double standardFrequency = std::stod(standardFrequencyString);
                table106[standardFrequency] = voltageList;
            }
            temp = jf["Table 107"];
            for (auto const& [standardFrequencyString, voltageList] : temp)
            {
                double standardFrequency = std::stod(standardFrequencyString);
                table107[standardFrequency] = voltageList;
            }
            temp = jf["Table 108"];
            for (auto const& [standardFrequencyString, voltageList] : temp)
            {
                double standardFrequency = std::stod(standardFrequencyString);
                table108[standardFrequency] = voltageList;
            }
            temp = jf["Table 109"];
            for (auto const& [standardFrequencyString, voltageList] : temp)
            {
                double standardFrequency = std::stod(standardFrequencyString);
                table109[standardFrequency] = voltageList;
            }
            temp = jf["Table 110"];
            for (auto const& [standardFrequencyString, voltageList] : temp)
            {
                double standardFrequency = std::stod(standardFrequencyString);
                table110[standardFrequency] = voltageList;
            }
        }
    }

    InsulationIEC61558Model(json data) {
        _data = data;
        table14 = data["IEC_61558-1"]["Table 14"];
        table20 = data["IEC_61558-1"]["Table 20"];
        table21 = data["IEC_61558-1"]["Table 21"];
        table22 = data["IEC_61558-1"]["Table 22"];
        table23 = data["IEC_61558-1"]["Table 23"];
        table102 = data["IEC_61558-2-16"]["Table 102"];
        table103 = data["IEC_61558-2-16"]["Table 103"];
        table104 = data["IEC_61558-2-16"]["Table 104"];
        std::map<std::string, std::vector<std::pair<double, double>>> temp;
        temp = data["IEC_61558-2-16"]["Table 105"];
        for (auto const& [standardFrequencyString, voltageList] : temp)
        {
            double standardFrequency = std::stod(standardFrequencyString);
            table105[standardFrequency] = voltageList;
        }
        temp = data["IEC_61558-2-16"]["Table 106"];
        for (auto const& [standardFrequencyString, voltageList] : temp)
        {
            double standardFrequency = std::stod(standardFrequencyString);
            table106[standardFrequency] = voltageList;
        }
        temp = data["IEC_61558-2-16"]["Table 107"];
        for (auto const& [standardFrequencyString, voltageList] : temp)
        {
            double standardFrequency = std::stod(standardFrequencyString);
            table107[standardFrequency] = voltageList;
        }
        temp = data["IEC_61558-2-16"]["Table 108"];
        for (auto const& [standardFrequencyString, voltageList] : temp)
        {
            double standardFrequency = std::stod(standardFrequencyString);
            table108[standardFrequency] = voltageList;
        }
        temp = data["IEC_61558-2-16"]["Table 109"];
        for (auto const& [standardFrequencyString, voltageList] : temp)
        {
            double standardFrequency = std::stod(standardFrequencyString);
            table109[standardFrequency] = voltageList;
        }
        temp = data["IEC_61558-2-16"]["Table 110"];
        for (auto const& [standardFrequencyString, voltageList] : temp)
        {
            double standardFrequency = std::stod(standardFrequencyString);
            table110[standardFrequency] = voltageList;
        }
    }
};

class InsulationIEC60335Model : public InsulationStandard {
  private:
    json _data;
  public:
    double calculate_distance_through_insulation(InputsWrapper inputs);
    double calculate_withstand_voltage(InputsWrapper inputs);
    double calculate_clearance(InputsWrapper inputs);
    double calculate_creepage_distance(InputsWrapper inputs, bool includeClearance = false);
    double get_rated_impulse_withstand_voltage(OvervoltageCategory overvoltageCategory, double ratedVoltage);
    double get_clearance_table_16(PollutionDegree pollutionDegree, WiringTechnology wiringType, InsulationType insulationType, double ratedImpulseWithstandVoltage);
    double get_distance_through_insulation_table_19(OvervoltageCategory overvoltageCategory, double ratedVoltage);
    double get_withstand_voltage_table_7(InsulationType insulationType, double ratedVoltage);
    double get_withstand_voltage_formula_table_7(InsulationType insulationType, double workingVoltage);
    double get_creepage_distance_table_17(Cti cti, PollutionDegree pollutionDegree, double workingVoltage);
    double get_creepage_distance_table_18(Cti cti, PollutionDegree pollutionDegree, double workingVoltage);

    double ieC60335MaximumStandardFrequency = 30000;
    std::map<std::string, std::vector<std::pair<double, double>>> table7;
    std::map<std::string, std::vector<std::pair<double, double>>> table15;
    std::vector<std::pair<double, double>> table16;
    std::map<std::string, std::map<std::string, std::vector<std::pair<double, double>>>> table17;
    std::map<std::string, std::map<std::string, std::vector<std::pair<double, double>>>> table18;
    std::map<std::string, std::vector<std::pair<double, double>>> table19;


    InsulationIEC60335Model() {
        std::string filePath = __FILE__;
        auto masPath = filePath.substr(0, filePath.rfind("/"));
        {
            auto dataFilePath = masPath + "/data/insulation_standards/IEC_60335-1.json";
            std::ifstream jsonFile(dataFilePath);
            json jf = json::parse(jsonFile);
            table7 = jf["Table 7"];
            table15 = jf["Table 15"];
            table16 = jf["Table 16"];
            table17 = jf["Table 17"];
            table18 = jf["Table 18"];
            table19 = jf["Table 19"];
        }
    }

    InsulationIEC60335Model(json data) {
        _data = data;
        table7 = data["IEC_60335-1"]["Table 7"];
        table15 = data["IEC_60335-1"]["Table 15"];
        table16 = data["IEC_60335-1"]["Table 16"];
        table17 = data["IEC_60335-1"]["Table 17"];
        table18 = data["IEC_60335-1"]["Table 18"];
        table19 = data["IEC_60335-1"]["Table 19"];
    }
};

class InsulationCoordinator {
  private:
  protected:
    std::shared_ptr<InsulationIEC60664Model> _insulationIEC60664Model;
    std::shared_ptr<InsulationIEC62368Model> _insulationIEC62368Model;
    std::shared_ptr<InsulationIEC61558Model> _insulationIEC61558Model;
    std::shared_ptr<InsulationIEC60335Model> _insulationIEC60335Model;

  public:
    double calculate_withstand_voltage(InputsWrapper inputs);
    double calculate_clearance(InputsWrapper inputs);
    double calculate_creepage_distance(InputsWrapper inputs, bool includeClearance = false);
    double calculate_distance_through_insulation(InputsWrapper inputs);
    InsulationCoordinationOutput calculate_insulation_coordination(InputsWrapper inputs);
    CoilSectionInterface calculate_coil_section_interface_layers(InputsWrapper inputs, WireWrapper leftWire, WireWrapper rightWire, InsulationMaterialWrapper insulationMaterial);
    static bool can_fully_insulated_wire_be_used(InputsWrapper inputs);

    InsulationCoordinator() {
        _insulationIEC60664Model = std::make_shared<InsulationIEC60664Model>(InsulationIEC60664Model());
        _insulationIEC62368Model = std::make_shared<InsulationIEC62368Model>(InsulationIEC62368Model());
        _insulationIEC61558Model = std::make_shared<InsulationIEC61558Model>(InsulationIEC61558Model());
        _insulationIEC60335Model = std::make_shared<InsulationIEC60335Model>(InsulationIEC60335Model());
    }
    InsulationCoordinator(json data) {
        _insulationIEC60664Model = std::make_shared<InsulationIEC60664Model>(InsulationIEC60664Model(data));
        _insulationIEC62368Model = std::make_shared<InsulationIEC62368Model>(InsulationIEC62368Model(data));
        _insulationIEC61558Model = std::make_shared<InsulationIEC61558Model>(InsulationIEC61558Model(data));
        _insulationIEC60335Model = std::make_shared<InsulationIEC60335Model>(InsulationIEC60335Model(data));
    }

    virtual ~InsulationCoordinator() = default;
};

} // namespace OpenMagnetics