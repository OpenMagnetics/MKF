#include "Insulation.h"
#include <UnitTest++.h>
#include "TestingUtils.h"
#include "Utils.h"


SUITE(Insulation) {
    // auto standardCoordinator = OpenMagnetics::InsulationCoordinator();
    // auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
    // auto cti = OpenMagnetics::Cti::GROUP_I;
    // double maximumVoltageRms = 666;
    // double maximumVoltagePeak = 800;
    // double frequency = 30000;
    // OpenMagnetics::DimensionWithTolerance altitude;
    // OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;

    TEST(IEC_60664_Load_Data) {
        auto standard = OpenMagnetics::InsulationIEC60664Model();
        OpenMagnetics::DimensionWithTolerance altitude;
        altitude.set_maximum(2000);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641};
        double maximumVoltageRms = 666;
        double maximumVoltagePeak = 800;
        double frequency = 30000;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0024, creepageDistance);
    }

    TEST(Test_Coordinated_Creepage_Distance) {
        double maximumVoltageRms = 666;
        double maximumVoltagePeak = 800;
        OpenMagnetics::DimensionWithTolerance altitude;
        OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        double frequency = 30000;
        auto standardCoordinator = OpenMagnetics::InsulationCoordinator();
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641, OpenMagnetics::InsulationStandards::IEC_623681};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standardCoordinator.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0024, creepageDistance);
    }

    TEST(Test_Coordinated_Clearance) {
        double maximumVoltageRms = 666;
        double maximumVoltagePeak = 800;
        auto cti = OpenMagnetics::Cti::GROUP_I;
        OpenMagnetics::DimensionWithTolerance altitude;
        OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;
        double frequency = 30000;
        auto standardCoordinator = OpenMagnetics::InsulationCoordinator();
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641, OpenMagnetics::InsulationStandards::IEC_623681};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standardCoordinator.calculate_clearance(inputs);
        CHECK_EQUAL(0.003, clearance);
    }

    TEST(Test_Insulation_Web_0) {
        auto standardCoordinator = OpenMagnetics::InsulationCoordinator();
        std::string inputString = R"({"designRequirements":{"insulation":{"altitude":{"maximum":2000},"cti":"Group I","pollutionDegree":"P1","overvoltageCategory":"OVC-I","insulationType":"Basic","mainSupplyVoltage":{"maximum":400},"standards":["IEC 60664-1"]},"magnetizingInductance":{"nominal":0.00001},"name":"My Design Requirements","turnsRatios":[],"wiringTechnology":"Wound"},"operatingPoints":[{"conditions":{"ambientTemperature":25},"excitationsPerWinding":[{"frequency":30000,"voltage":{"processed":{"dutyCycle":0.5,"peak":800,"peakToPeak":1600,"rms":666,"offset":0,"label":"Rectangular"}}}]}]})";
        OpenMagnetics::InputsWrapper inputs = OpenMagnetics::InputsWrapper(json::parse(inputString), false);

        auto clearance = standardCoordinator.calculate_clearance(inputs);
        auto creepageDistance = standardCoordinator.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.003, clearance);
        CHECK_EQUAL(0.0024, creepageDistance);
    }

    TEST(Test_Insulation_Web_1) {
        auto standardCoordinator = OpenMagnetics::InsulationCoordinator();
        std::string inputString = R"({"designRequirements":{"insulation":{"altitude":{"maximum":2000},"cti":"Group I","pollutionDegree":"P1","overvoltageCategory":"OVC-I","insulationType":"Basic","mainSupplyVoltage":{"maximum":400},"standards":["IEC 60664-1"]},"magnetizingInductance":{"nominal":0.00001},"name":"My Design Requirements","turnsRatios":[],"wiringTechnology":"Wound"},"operatingPoints":[{"conditions":{"ambientTemperature":25},"excitationsPerWinding":[{"frequency":100000,"voltage":{"processed":{"dutyCycle":0.5,"peak":800,"peakToPeak":1600,"rms":666,"offset":0,"label":"Rectangular"}}}]}]})";
        OpenMagnetics::InputsWrapper inputs = OpenMagnetics::InputsWrapper(json::parse(inputString), false);

        auto clearance = standardCoordinator.calculate_clearance(inputs);
        auto creepageDistance = standardCoordinator.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.003, clearance);
        CHECK_EQUAL(0.0024, creepageDistance);
    }

    TEST(Test_Insulation_Web_2) {
        auto standardCoordinator = OpenMagnetics::InsulationCoordinator();
        std::string inputString = R"({"designRequirements":{"insulation":{"altitude":{"excludeMaximum":null,"excludeMinimum":null,"maximum":2000.0,"minimum":null,"nominal":null},"cti":"Group IIIB","insulationType":"Double","mainSupplyVoltage":{"excludeMaximum":null,"excludeMinimum":null,"maximum":400.0,"minimum":null,"nominal":null},"overvoltageCategory":"OVC-III","pollutionDegree":"P2","standards":["IEC 60664-1"]},"isolationSides":null,"leakageInductance":null,"magnetizingInductance":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0001},"market":null,"maximumDimensions":null,"maximumWeight":null,"name":"My Design Requirements","operatingTemperature":null,"strayCapacitance":null,"terminalType":null,"topology":null,"turnsRatios":[{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":1.0}],"wiringTechnology":null},"operatingPoints":[{"conditions":{"ambientRelativeHumidity":null,"ambientTemperature":42.0,"cooling":null,"name":null},"excitationsPerWinding":[{"current":{"harmonics":{"amplitudes":[1.1608769501236793e-14,4.05366124583194,1.787369544444173e-15,0.4511310569983995,9.749053004706756e-16,0.16293015292554872,4.036157626725542e-16,0.08352979924600704,3.4998295008010614e-16,0.0508569581336163,3.1489164048780735e-16,0.034320410449418075,3.142469873118059e-16,0.024811988673843106,2.3653352035940994e-16,0.018849001010678823,2.9306524147249266e-16,0.014866633059596499,1.796485796132569e-16,0.012077180559557785,1.6247782523152451e-16,0.010049063750920609,1.5324769149805092e-16,0.008529750975091871,1.0558579733068502e-16,0.007363501410705499,7.513269775674661e-17,0.006450045785294609,5.871414177162291e-17,0.005722473794997712,9.294731722001391e-17,0.005134763398167541,1.194820309200107e-16,0.004654430423785411,8.2422739080512e-17,0.004258029771397032,9.5067306351894e-17,0.0039283108282380024,1.7540347128474968e-16,0.0036523670873925395,9.623794010508822e-17,0.0034204021424253787,1.4083470894369491e-16,0.0032248884817922927,1.4749333016985644e-16,0.0030599828465501895,1.0448590642474364e-16,0.002921112944200383,7.575487373767413e-18,0.002804680975178716,7.419510610361002e-17,0.0027078483284668376,3.924741709073613e-17,0.0026283777262804953,2.2684279102637236e-17,0.0025645167846443107,8.997077625295079e-17,0.0025149120164513483,7.131074184849219e-17,0.0024785457043284276,9.354417496250849e-17,0.0024546904085875065,1.2488589642405877e-16,0.0024428775264784264],"frequencies":[0.0,100000.0,200000.0,300000.0,400000.0,500000.0,600000.0,700000.0,800000.0,900000.0,1000000.0,1100000.0,1200000.0,1300000.0,1400000.0,1500000.0,1600000.0,1700000.0,1800000.0,1900000.0,2000000.0,2100000.0,2200000.0,2300000.0,2400000.0,2500000.0,2600000.0,2700000.0,2800000.0,2900000.0,3000000.0,3100000.0,3200000.0,3300000.0,3400000.0,3500000.0,3600000.0,3700000.0,3800000.0,3900000.0,4000000.0,4100000.0,4200000.0,4300000.0,4400000.0,4500000.0,4600000.0,4700000.0,4800000.0,4900000.0,5000000.0,5100000.0,5200000.0,5300000.0,5400000.0,5500000.0,5600000.0,5700000.0,5800000.0,5900000.0,6000000.0,6100000.0,6200000.0,6300000.0]},"processed":{"acEffectiveFrequency":110746.40291779551,"average":null,"dutyCycle":0.5,"effectiveFrequency":110746.40291779551,"label":"Triangular","offset":0.0,"peak":5.0,"peakToPeak":10.0,"phase":null,"rms":2.8874560332150576,"thd":0.12151487440704967},"waveform":{"ancillaryLabel":null,"data":[-5.0,5.0,-5.0],"numberPeriods":null,"time":[0.0,5e-06,1e-05]}},"frequency":100000.0,"magneticFieldStrength":null,"magneticFluxDensity":null,"magnetizingCurrent":null,"name":"Primary winding excitation","voltage":{"harmonics":{"amplitudes":[24.2890625,57.92076613061847,1.421875,19.27588907896988,1.421875,11.528257939603122,1.421875,8.194467538528329,1.421875,6.331896912839248,1.421875,5.137996046859012,1.421875,4.304077056139349,1.421875,3.6860723299088454,1.421875,3.207698601961777,1.421875,2.8247804703632298,1.421875,2.509960393415185,1.421875,2.2453859950684323,1.421875,2.01890737840567,1.421875,1.8219644341144872,1.421875,1.6483482744897402,1.421875,1.4934420157473332,1.421875,1.3537375367153817,1.421875,1.2265178099275544,1.421875,1.1096421410704556,1.421875,1.0013973584174929,1.421875,0.9003924136274832,1.421875,0.8054822382572133,1.421875,0.7157117294021269,1.421875,0.6302738400635857,1.421875,0.5484777114167545,1.421875,0.46972405216147894,1.421875,0.3934858059809043,1.421875,0.31929270856030145,1.421875,0.24671871675852053,1.421875,0.17537155450693565,1.421875,0.10488380107099537,1.421875,0.034905072061178544],"frequencies":[0.0,100000.0,200000.0,300000.0,400000.0,500000.0,600000.0,700000.0,800000.0,900000.0,1000000.0,1100000.0,1200000.0,1300000.0,1400000.0,1500000.0,1600000.0,1700000.0,1800000.0,1900000.0,2000000.0,2100000.0,2200000.0,2300000.0,2400000.0,2500000.0,2600000.0,2700000.0,2800000.0,2900000.0,3000000.0,3100000.0,3200000.0,3300000.0,3400000.0,3500000.0,3600000.0,3700000.0,3800000.0,3900000.0,4000000.0,4100000.0,4200000.0,4300000.0,4400000.0,4500000.0,4600000.0,4700000.0,4800000.0,4900000.0,5000000.0,5100000.0,5200000.0,5300000.0,5400000.0,5500000.0,5600000.0,5700000.0,5800000.0,5900000.0,6000000.0,6100000.0,6200000.0,6300000.0]},"processed":{"acEffectiveFrequency":591485.5360118389,"average":null,"dutyCycle":0.5,"effectiveFrequency":553357.3374711702,"label":"Rectangular","offset":0.0,"peak":70.5,"peakToPeak":100.0,"phase":null,"rms":51.572309672924284,"thd":0.4833151484524849},"waveform":{"ancillaryLabel":null,"data":[-20.5,70.5,70.5,-20.5,-20.5],"numberPeriods":null,"time":[0.0,0.0,5e-06,5e-06,1e-05]}}},{"current":{"harmonics":null,"processed":{"acEffectiveFrequency":null,"average":null,"dutyCycle":0.5,"effectiveFrequency":null,"label":"Triangular","offset":0.0,"peak":null,"peakToPeak":10.0,"phase":null,"rms":null,"thd":null},"waveform":{"ancillaryLabel":null,"data":[-5.0,5.0,-5.0],"numberPeriods":null,"time":[0.0,5e-06,1e-05]}},"frequency":100000.0,"magneticFieldStrength":null,"magneticFluxDensity":null,"magnetizingCurrent":null,"name":"Primary winding excitation","voltage":{"harmonics":null,"processed":{"acEffectiveFrequency":null,"average":null,"dutyCycle":0.5,"effectiveFrequency":null,"label":"Rectangular","offset":0.0,"peak":null,"peakToPeak":100.0,"phase":null,"rms":null,"thd":null},"waveform":{"ancillaryLabel":null,"data":[-20.5,70.5,70.5,-20.5,-20.5],"numberPeriods":null,"time":[0.0,0.0,5e-06,5e-06,1e-05]}}}],"name":"Operating Point No. 1"}]})";
        OpenMagnetics::InputsWrapper inputs = OpenMagnetics::InputsWrapper(json::parse(inputString), false);

        auto clearance = standardCoordinator.calculate_clearance(inputs);
    }
}

SUITE(CoilSectionsInterface) {
    // auto standardCoordinator = OpenMagnetics::InsulationCoordinator();
    // auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
    // auto cti = OpenMagnetics::Cti::GROUP_I;
    // double maximumVoltageRms = 1000;
    // double maximumVoltagePeak = 1800;
    // double frequency = 30000;
    // OpenMagnetics::DimensionWithTolerance altitude;
    // OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;
    // auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;

    TEST(Test_Basic_SIW_SIW_OVC_I_Kapton) {
        auto standardCoordinator = OpenMagnetics::InsulationCoordinator();
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto cti = OpenMagnetics::Cti::GROUP_I;
        double maximumVoltageRms = 1000;
        double maximumVoltagePeak = 1800;
        double frequency = 30000;
        OpenMagnetics::DimensionWithTolerance altitude;
        OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(800);
        overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(1);
        inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
        inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<OpenMagnetics::IsolationSide>{OpenMagnetics::IsolationSide::PRIMARY, OpenMagnetics::IsolationSide::PRIMARY});
        auto insulationMaterial = OpenMagnetics::find_insulation_material_by_name("Kapton HN");
        auto leftWire = OpenMagnetics::find_wire_by_name("Litz SXXL825/44FX-3(MWXX)");
        auto rightWire = OpenMagnetics::find_wire_by_name("Litz SXXL825/44FX-3(MWXX)");

        auto coilSectionInterface = standardCoordinator.calculate_coil_section_interface_layers(inputs, leftWire,  rightWire, insulationMaterial).value();
        CHECK(coilSectionInterface.get_total_margin_tape_distance() == 0);
        CHECK_EQUAL(1UL, coilSectionInterface.get_number_layers_insulation());
        CHECK(OpenMagnetics::CoilSectionInterface::LayerPurpose::INSULATING == coilSectionInterface.get_layer_purpose());
    }

    TEST(Test_Reinforced_SIW_SIW_OVC_I_Tecroll_10B) {
        auto standardCoordinator = OpenMagnetics::InsulationCoordinator();
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto cti = OpenMagnetics::Cti::GROUP_I;
        double maximumVoltageRms = 1000;
        double maximumVoltagePeak = 1800;
        double frequency = 30000;
        OpenMagnetics::DimensionWithTolerance altitude;
        OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(800);
        overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(1);
        inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
        inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<OpenMagnetics::IsolationSide>{OpenMagnetics::IsolationSide::PRIMARY, OpenMagnetics::IsolationSide::PRIMARY});
        auto insulationMaterial = OpenMagnetics::find_insulation_material_by_name("Tecroll 10B");
        auto leftWire = OpenMagnetics::find_wire_by_name("Litz SXXL825/44FX-3(MWXX)");
        auto rightWire = OpenMagnetics::find_wire_by_name("Litz SXXL825/44FX-3(MWXX)");

        auto coilSectionInterface = standardCoordinator.calculate_coil_section_interface_layers(inputs, leftWire,  rightWire, insulationMaterial).value();
        CHECK(coilSectionInterface.get_total_margin_tape_distance() > 0);
        CHECK_EQUAL(2UL, coilSectionInterface.get_number_layers_insulation());
        CHECK(OpenMagnetics::CoilSectionInterface::LayerPurpose::INSULATING == coilSectionInterface.get_layer_purpose());
    }

    TEST(Test_Reinforced_SIW_SIW_OVC_IV_Kapton) {
        auto standardCoordinator = OpenMagnetics::InsulationCoordinator();
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto cti = OpenMagnetics::Cti::GROUP_I;
        double maximumVoltageRms = 1000;
        double maximumVoltagePeak = 1800;
        double frequency = 30000;
        OpenMagnetics::DimensionWithTolerance altitude;
        OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(800);
        overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(1);
        inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
        inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<OpenMagnetics::IsolationSide>{OpenMagnetics::IsolationSide::PRIMARY, OpenMagnetics::IsolationSide::PRIMARY});
        auto insulationMaterial = OpenMagnetics::find_insulation_material_by_name("Kapton HN");
        auto leftWire = OpenMagnetics::find_wire_by_name("Litz SXXL825/44FX-3(MWXX)");
        auto rightWire = OpenMagnetics::find_wire_by_name("Litz SXXL825/44FX-3(MWXX)");

        auto coilSectionInterface = standardCoordinator.calculate_coil_section_interface_layers(inputs, leftWire,  rightWire, insulationMaterial).value();
        CHECK(coilSectionInterface.get_total_margin_tape_distance() > 0);
        CHECK_EQUAL(1UL, coilSectionInterface.get_number_layers_insulation());
        CHECK(OpenMagnetics::CoilSectionInterface::LayerPurpose::INSULATING == coilSectionInterface.get_layer_purpose());
    }

    TEST(Test_Basic_SIW_SIW_OVC_IV_Kapton) {
        auto standardCoordinator = OpenMagnetics::InsulationCoordinator();
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto cti = OpenMagnetics::Cti::GROUP_I;
        double maximumVoltageRms = 1000;
        double maximumVoltagePeak = 1800;
        double frequency = 30000;
        OpenMagnetics::DimensionWithTolerance altitude;
        OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(800);
        overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(1);
        inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
        inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<OpenMagnetics::IsolationSide>{OpenMagnetics::IsolationSide::PRIMARY, OpenMagnetics::IsolationSide::PRIMARY});
        auto insulationMaterial = OpenMagnetics::find_insulation_material_by_name("Kapton HN");
        auto leftWire = OpenMagnetics::find_wire_by_name("Litz SXXL825/44FX-3(MWXX)");
        auto rightWire = OpenMagnetics::find_wire_by_name("Litz SXXL825/44FX-3(MWXX)");

        auto coilSectionInterface = standardCoordinator.calculate_coil_section_interface_layers(inputs, leftWire,  rightWire, insulationMaterial).value();
        CHECK(coilSectionInterface.get_total_margin_tape_distance() > 0);
        CHECK_EQUAL(1UL, coilSectionInterface.get_number_layers_insulation());
        CHECK(OpenMagnetics::CoilSectionInterface::LayerPurpose::INSULATING == coilSectionInterface.get_layer_purpose());
    }

    TEST(Test_Basic_SIW_SIW_OVC_IV_ETFE) {
        auto standardCoordinator = OpenMagnetics::InsulationCoordinator();
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto cti = OpenMagnetics::Cti::GROUP_I;
        double maximumVoltageRms = 1000;
        double maximumVoltagePeak = 1800;
        double frequency = 30000;
        OpenMagnetics::DimensionWithTolerance altitude;
        OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(800);
        overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(1);
        inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
        inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<OpenMagnetics::IsolationSide>{OpenMagnetics::IsolationSide::PRIMARY, OpenMagnetics::IsolationSide::PRIMARY});
        auto insulationMaterial = OpenMagnetics::find_insulation_material_by_name("ETFE");
        auto leftWire = OpenMagnetics::find_wire_by_name("Litz SXXL825/44FX-3(MWXX)");
        auto rightWire = OpenMagnetics::find_wire_by_name("Litz SXXL825/44FX-3(MWXX)");

        auto coilSectionInterface = standardCoordinator.calculate_coil_section_interface_layers(inputs, leftWire,  rightWire, insulationMaterial).value();
        CHECK(coilSectionInterface.get_total_margin_tape_distance() > 0);
        CHECK_EQUAL(3UL, coilSectionInterface.get_number_layers_insulation());
        CHECK(OpenMagnetics::CoilSectionInterface::LayerPurpose::INSULATING == coilSectionInterface.get_layer_purpose());
    }

    TEST(Test_Basic_DIW_SIW_OVC_I_ETFE) {
        auto standardCoordinator = OpenMagnetics::InsulationCoordinator();
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto cti = OpenMagnetics::Cti::GROUP_I;
        double maximumVoltageRms = 1000;
        double maximumVoltagePeak = 1800;
        double frequency = 30000;
        OpenMagnetics::DimensionWithTolerance altitude;
        OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(800);
        overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(1);
        inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
        inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<OpenMagnetics::IsolationSide>{OpenMagnetics::IsolationSide::PRIMARY, OpenMagnetics::IsolationSide::PRIMARY});
        auto insulationMaterial = OpenMagnetics::find_insulation_material_by_name("ETFE");
        auto leftWire = OpenMagnetics::find_wire_by_name("Litz SXXL20/34FX-3(MWXX)");
        auto rightWire = OpenMagnetics::find_wire_by_name("Litz DXXL07/28TXX-3(MWXX)");

        auto coilSectionInterface = standardCoordinator.calculate_coil_section_interface_layers(inputs, leftWire,  rightWire, insulationMaterial).value();
        CHECK(coilSectionInterface.get_total_margin_tape_distance() == 0);
        CHECK_EQUAL(1UL, coilSectionInterface.get_number_layers_insulation());
        CHECK(OpenMagnetics::CoilSectionInterface::LayerPurpose::MECHANICAL == coilSectionInterface.get_layer_purpose());
    }

    TEST(Test_Basic_DIW_Enammeled_Wire_OVC_I_ETFE) {
        auto standardCoordinator = OpenMagnetics::InsulationCoordinator();
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        OpenMagnetics::DimensionWithTolerance altitude;
        OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;
        standardCoordinator = OpenMagnetics::InsulationCoordinator();
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(800);
        overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto cti = OpenMagnetics::Cti::GROUP_I;
        double maximumVoltageRms = 1000;
        double maximumVoltagePeak = 1800;
        double frequency = 30000;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(1);
        inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
        inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<OpenMagnetics::IsolationSide>{OpenMagnetics::IsolationSide::PRIMARY, OpenMagnetics::IsolationSide::PRIMARY});
        auto insulationMaterial = OpenMagnetics::find_insulation_material_by_name("ETFE");
        auto leftWire = OpenMagnetics::find_wire_by_name("Round 0.016 - Grade 1");
        auto rightWire = OpenMagnetics::find_wire_by_name("Litz DXXL07/28TXX-3(MWXX)");

        auto coilSectionInterface = standardCoordinator.calculate_coil_section_interface_layers(inputs, leftWire,  rightWire, insulationMaterial).value();
        CHECK(coilSectionInterface.get_total_margin_tape_distance() == 0);
        CHECK_EQUAL(1UL, coilSectionInterface.get_number_layers_insulation());
        CHECK(OpenMagnetics::CoilSectionInterface::LayerPurpose::INSULATING == coilSectionInterface.get_layer_purpose());
    }

    TEST(Test_Basic_TIW_OVC_I_ETFE) {
        auto standardCoordinator = OpenMagnetics::InsulationCoordinator();
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        double maximumVoltageRms = 1000;
        double maximumVoltagePeak = 1800;
        double frequency = 30000;
        OpenMagnetics::DimensionWithTolerance altitude;
        OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(800);
        overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(1);
        inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
        inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<OpenMagnetics::IsolationSide>{OpenMagnetics::IsolationSide::PRIMARY, OpenMagnetics::IsolationSide::PRIMARY});
        auto insulationMaterial = OpenMagnetics::find_insulation_material_by_name("ETFE");
        auto leftWire = OpenMagnetics::find_wire_by_name("Round T28A01TXXX-1.5");
        auto rightWire = OpenMagnetics::find_wire_by_name("Round 0.016 - Grade 1");

        auto coilSectionInterface = standardCoordinator.calculate_coil_section_interface_layers(inputs, leftWire,  rightWire, insulationMaterial).value();
        CHECK(coilSectionInterface.get_total_margin_tape_distance() == 0);
        CHECK_EQUAL(1UL, coilSectionInterface.get_number_layers_insulation());
        CHECK(OpenMagnetics::CoilSectionInterface::LayerPurpose::MECHANICAL == coilSectionInterface.get_layer_purpose());
    }

    TEST(Test_Basic_FIW_OVC_I_ETFE) {
        auto standardCoordinator = OpenMagnetics::InsulationCoordinator();
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto cti = OpenMagnetics::Cti::GROUP_I;
        double maximumVoltageRms = 1000;
        double maximumVoltagePeak = 1800;
        double frequency = 30000;
        OpenMagnetics::DimensionWithTolerance altitude;
        OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_623681};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(1);
        inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
        inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<OpenMagnetics::IsolationSide>{OpenMagnetics::IsolationSide::PRIMARY, OpenMagnetics::IsolationSide::PRIMARY});
        auto insulationMaterial = OpenMagnetics::find_insulation_material_by_name("ETFE");
        auto leftWire = OpenMagnetics::find_wire_by_name("Round 0.071 - FIW 9");
        auto rightWire = OpenMagnetics::find_wire_by_name("Round 0.016 - Grade 1");

        auto coilSectionInterface = standardCoordinator.calculate_coil_section_interface_layers(inputs, leftWire,  rightWire, insulationMaterial).value();
        CHECK(coilSectionInterface.get_total_margin_tape_distance() == 0);
        CHECK_EQUAL(1UL, coilSectionInterface.get_number_layers_insulation());
        CHECK(OpenMagnetics::CoilSectionInterface::LayerPurpose::MECHANICAL == coilSectionInterface.get_layer_purpose());
    }
}

SUITE(CreepageDistance_IEC_60664) {

    auto standard = OpenMagnetics::InsulationIEC60664Model();
    auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641};
    auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
    double maximumVoltageRms = 666;
    double maximumVoltagePeak = 800;
    double frequency = 30000;
    OpenMagnetics::DimensionWithTolerance altitude;
    OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;


    TEST(Creepage_Distance_Basic_P1_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0024, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0048, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.004, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.008, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.01, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.02, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0024, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0048, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0056, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0112, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.011, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.022, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0024, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0048, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.008, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.016, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0125, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.025, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_I_High_Frequency) {
        frequency = 700000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0038, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_I_High_Frequency) {
        frequency = 700000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0076, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_I_High_Frequency) {
        frequency = 700000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.00456, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_I_High_Frequency) {
        frequency = 700000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.00912, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_I_High_Frequency) {
        frequency = 700000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.01, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_I_High_Frequency) {
        frequency = 700000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.02, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_II_High_Frequency) {
        frequency = 700000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0038, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_II_High_Frequency) {
        frequency = 700000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0076, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_II_High_Frequency) {
        frequency = 700000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0056, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_II_High_Frequency) {
        frequency = 700000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0112, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_II_High_Frequency) {
        frequency = 700000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.011, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_II_High_Frequency) {
        frequency = 700000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.022, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_IIIA_High_Frequency) {
        frequency = 700000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0038, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_IIIA_High_Frequency) {
        frequency = 700000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0076, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_IIIA_High_Frequency) {
        frequency = 700000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.008, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_IIIA_High_Frequency) {
        frequency = 700000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.016, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_IIIA_High_Frequency) {
        frequency = 700000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0125, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_IIIA_High_Frequency) {
        frequency = 700000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.025, creepageDistance);
    }
}

SUITE(Clearance_IEC_60664) {
    auto standard = OpenMagnetics::InsulationIEC60664Model();
    auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641};
    auto cti = OpenMagnetics::Cti::GROUP_I;
    double maximumVoltageRms = 69;
    double maximumVoltagePeak = 260;
    double frequency = 30000;
    OpenMagnetics::DimensionWithTolerance altitude;
    OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;

    TEST(Clearance_Basic_P1_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.00004, clearance);
    }

    TEST(Clearance_Reinforced_P1_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0001, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0002, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0002, clearance);
    }

    TEST(Clearance_Basic_P3_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0008, clearance);
    }

    TEST(Clearance_Reinforced_P3_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0008, clearance);
    }

    TEST(Clearance_Basic_P1_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0001, clearance);
    }

    TEST(Clearance_Reinforced_P1_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0005, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0002, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0005, clearance);
    }

    TEST(Clearance_Basic_P3_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0008, clearance);
    }

    TEST(Clearance_Reinforced_P3_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0008, clearance);
    }

    TEST(Clearance_Basic_P1_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0005, clearance);
    }

    TEST(Clearance_Reinforced_P1_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0015, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0005, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0015, clearance);
    }

    TEST(Clearance_Basic_P3_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0008, clearance);
    }

    TEST(Clearance_Reinforced_P3_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0015, clearance);
    }

    TEST(Clearance_Basic_P1_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0015, clearance);
    }

    TEST(Clearance_Reinforced_P1_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.003, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0015, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.003, clearance);
    }

    TEST(Clearance_Basic_P3_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0015, clearance);
    }

    TEST(Clearance_Reinforced_P3_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.003, clearance);
    }

    TEST(Clearance_Basic_P1_OVC_I_High_Altitude_Low_Frequency) {
        altitude.set_maximum(8000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.00004 * 2.25, clearance);
    }

    TEST(Clearance_Basic_P1_OVC_I_Low_Altitude_Low_Frequency_High_Voltage) {
        maximumVoltageRms = 666;
        maximumVoltagePeak = 800;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.003, clearance);
    }

    TEST(Clearance_Basic_P1_OVC_I_Low_Altitude_High_Frequency_High_Voltage) {
        frequency = 500000;
        maximumVoltageRms = 666;
        maximumVoltagePeak = 800;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.003, clearance);
    }
}

SUITE(Distance_Through_Insulation_IEC_60664) {
    auto standard = OpenMagnetics::InsulationIEC60664Model();
    auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641};
    auto cti = OpenMagnetics::Cti::GROUP_I;
    double maximumVoltageRms = 666;
    double maximumVoltagePeak = 800;
    double frequency = 30000;
    OpenMagnetics::DimensionWithTolerance altitude;
    OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;

    TEST(Distance_Through_Insulation_High_Frequency) {
        frequency = 500000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto dti = standard.calculate_distance_through_insulation(inputs);
        CHECK_EQUAL(0.00025, dti);
    }
}

SUITE(CreepageDistance_IEC_62368) {

    auto standard = OpenMagnetics::InsulationIEC62368Model();
    auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_623681};
    auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
    double maximumVoltageRms = 666;
    double maximumVoltagePeak = 800;
    double frequency = 30000;
    OpenMagnetics::DimensionWithTolerance altitude;
    OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;


    TEST(Creepage_Distance_Basic_P1_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.002, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0039, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0034, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0068, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0085, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0169, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.002, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0039, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0048, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0095, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0095, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0189, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.002, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0039, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0067, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0134, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0106, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0211, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_I_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.002, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_I_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0039, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_I_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0034, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_I_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0068, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_I_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0085, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_I_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0169, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_II_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.002, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_II_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0039, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_II_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0048, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_II_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0095, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_II_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0095, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_II_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0189, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_IIIA_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.002, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_IIIA_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0039, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_IIIA_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0067, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_IIIA_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0134, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_IIIA_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0106, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_IIIA_High_Frequency) {
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0211, creepageDistance);
    }
}

SUITE(Clearance_IEC_62368) {
    auto standard = OpenMagnetics::InsulationIEC62368Model();
    auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_623681};
    auto cti = OpenMagnetics::Cti::GROUP_I;
    double maximumVoltageRms = 666;
    double maximumVoltagePeak = 800;
    double frequency = 30000;
    OpenMagnetics::DimensionWithTolerance altitude;
    OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;

    TEST(Clearance_Basic_P1_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0018, clearance);
    }

    TEST(Clearance_Reinforced_P1_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0036, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0018, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0036, clearance);
    }

    TEST(Clearance_Basic_P3_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0018, clearance);
    }

    TEST(Clearance_Reinforced_P3_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0036, clearance);
    }

    TEST(Clearance_Basic_P1_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0018, clearance);
    }

    TEST(Clearance_Reinforced_P1_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0036, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0018, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0036, clearance);
    }

    TEST(Clearance_Basic_P3_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0018, clearance);
    }

    TEST(Clearance_Reinforced_P3_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0036, clearance);
    }

    TEST(Clearance_Basic_P1_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.003, clearance);
    }

    TEST(Clearance_Reinforced_P1_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0055, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.003, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0055, clearance);
    }

    TEST(Clearance_Basic_P3_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.003, clearance);
    }

    TEST(Clearance_Reinforced_P3_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0055, clearance);
    }

    TEST(Clearance_Basic_P1_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0055, clearance);
    }

    TEST(Clearance_Reinforced_P1_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.008, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0055, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.008, clearance);
    }

    TEST(Clearance_Basic_P3_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0055, clearance);
    }

    TEST(Clearance_Reinforced_P3_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.008, clearance);
    }

    TEST(Clearance_Basic_P1_OVC_I_High_Altitude_Low_Frequency) {
        altitude.set_maximum(5000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.00267, clearance);
    }

    TEST(Clearance_Reinforced_P1_OVC_I_High_Altitude_Low_Frequency) {
        altitude.set_maximum(5000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.00533, clearance);
    }    

    TEST(Clearance_Basic_P1_OVC_I_Low_Altitude_High_Frequency_High_Voltage_Peak) {
        altitude.set_maximum(2000);
        frequency = 400000;
        maximumVoltagePeak = 2000;
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0106, clearance);
    }

    TEST(Clearance_Reinforced_P1_OVC_I_Low_Altitude_High_Frequency_High_Voltage_Peak) {
        altitude.set_maximum(2000);
        frequency = 400000;
        maximumVoltagePeak = 2000;
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0212, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_I_Low_Altitude_High_Frequency_High_Voltage_Peak) {
        altitude.set_maximum(2000);
        frequency = 400000;
        maximumVoltagePeak = 2000;
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0132, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_I_Low_Altitude_High_Frequency_High_Voltage_Peak) {
        altitude.set_maximum(2000);
        frequency = 400000;
        maximumVoltagePeak = 2000;
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0264, clearance);
    }

    TEST(Clearance_Printed_Basic) {
        altitude.set_maximum(2000);
        frequency = 100000;
        maximumVoltagePeak = 2000;
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::PRINTED);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0038, clearance);
    }

    TEST(Clearance_Printed_Reinforced) {
        altitude.set_maximum(2000);
        frequency = 100000;
        maximumVoltagePeak = 2000;
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::PRINTED);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0044, clearance);
    }
}

SUITE(CreepageDistance_IEC_61558) {

    auto standard = OpenMagnetics::InsulationIEC61558Model();
    auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_615581};
    auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
    double maximumVoltageRms = 666;
    double maximumVoltagePeak = 800;
    double frequency = 30000;
    OpenMagnetics::DimensionWithTolerance altitude;
    OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;


    TEST(Creepage_Distance_Basic_P1_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.00195, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.00458, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.00342, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.00666, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0085, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.01749, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.00195, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.00458, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.00477, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0095, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0095, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.01899, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.00195, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.00458, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.00666, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.01332, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.01058, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.02132, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_High_Frequency) {
        frequency = 600000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.00290, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_High_Frequency) {
        frequency = 600000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0058, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_High_Frequency) {
        frequency = 600000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.00348, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_High_Frequency) {
        frequency = 600000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.00696, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_High_Frequency) {
        frequency = 600000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0085, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_High_Frequency) {
        frequency = 600000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.01749, creepageDistance);
    }
}

SUITE(Clearance_IEC_61558) {
    auto standard = OpenMagnetics::InsulationIEC61558Model();
    auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_615581};
    auto cti = OpenMagnetics::Cti::GROUP_I;
    double maximumVoltageRms = 666;
    double maximumVoltagePeak = 800;
    double frequency = 30000;
    OpenMagnetics::DimensionWithTolerance altitude;
    OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;

    TEST(Clearance_Basic_P1_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0, clearance);
    }

    TEST(Clearance_Reinforced_P1_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.003, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0055, clearance);
    }

    TEST(Clearance_Basic_P3_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.003, clearance);
    }

    TEST(Clearance_Reinforced_P3_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0055, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0055, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.008, clearance);
    }

    TEST(Clearance_Basic_P3_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0055, clearance);
    }

    TEST(Clearance_Reinforced_P3_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.008, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.008, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.014, clearance);
    }

    TEST(Clearance_Basic_P3_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.008, clearance);
    }

    TEST(Clearance_Reinforced_P3_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.014, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.014, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.025, clearance);
    }

    TEST(Clearance_Basic_P3_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.014, clearance);
    }

    TEST(Clearance_Reinforced_P3_OVC_IV_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.025, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_I_High_Altitude_Low_Frequency) {
        altitude.set_maximum(5000);
        mainSupplyVoltage.set_nominal(400);
        frequency = 30000;
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.00444, clearance);
    }

    TEST(Clearance_Printed_Basic) {
        altitude.set_maximum(2000);
        frequency = 30000;
        maximumVoltagePeak = 2000;
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::PRINTED);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0, clearance);
    }

    TEST(Clearance_Printed_Reinforced) {
        altitude.set_maximum(2000);
        frequency = 30000;
        maximumVoltagePeak = 2000;
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::PRINTED);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0, clearance);
    }
}

SUITE(Distance_Through_Insulation_IEC_61558) {
    auto standard = OpenMagnetics::InsulationIEC61558Model();
    auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_615581};
    auto cti = OpenMagnetics::Cti::GROUP_I;
    double maximumVoltageRms = 666;
    double maximumVoltagePeak = 800;
    double frequency = 30000;
    OpenMagnetics::DimensionWithTolerance altitude;
    OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;

    TEST(Distance_Through_Insulation_Basic_Solid) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto dti = standard.calculate_distance_through_insulation(inputs, false);
        CHECK_EQUAL(0, dti);
    }

    TEST(Distance_Through_Insulation_Supplementary_Solid) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::SUPPLEMENTARY;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto dti = standard.calculate_distance_through_insulation(inputs, false);
        CHECK_EQUAL(0.0008, dti);
    }

    TEST(Distance_Through_Insulation_Reinforced_Solid) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto dti = standard.calculate_distance_through_insulation(inputs, false);
        CHECK_EQUAL(0.00159, dti);
    }

    TEST(Distance_Through_Insulation_Basic_Thin_Sheet) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto dti = standard.calculate_distance_through_insulation(inputs);
        CHECK_EQUAL(0, dti);
    }

    TEST(Distance_Through_Insulation_Supplementary_Thin_Sheet) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::SUPPLEMENTARY;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto dti = standard.calculate_distance_through_insulation(inputs);
        CHECK_EQUAL(0.00021, dti);
    }

    TEST(Distance_Through_Insulation_Reinforced_Thin_Sheet) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto dti = standard.calculate_distance_through_insulation(inputs);
        CHECK_EQUAL(0.00042, dti);
    }

    TEST(Distance_Through_Insulation_High_Frequency) {
        frequency = 500000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto dti = standard.calculate_distance_through_insulation(inputs);
        CHECK_EQUAL(0.00025, dti);
    }
}

SUITE(Clearance_IEC_60335) {
    auto standard = OpenMagnetics::InsulationIEC60335Model();
    auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_603351};
    auto cti = OpenMagnetics::Cti::GROUP_I;
    double maximumVoltageRms = 250;
    double maximumVoltagePeak = 400;
    double frequency = 30000;
    OpenMagnetics::DimensionWithTolerance altitude;
    OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;

    TEST(Clearance_Basic_P1_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0005, clearance);
    }

    TEST(Clearance_Reinforced_P1_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0015, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0005, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0015, clearance);
    }

    TEST(Clearance_Basic_P3_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0008, clearance);
    }

    TEST(Clearance_Reinforced_P3_OVC_I_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0015, clearance);
    }

    TEST(Clearance_Basic_P1_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0015, clearance);
    }

    TEST(Clearance_Reinforced_P1_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.003, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0015, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.003, clearance);
    }

    TEST(Clearance_Basic_P3_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0015, clearance);
    }

    TEST(Clearance_Reinforced_P3_OVC_II_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.003, clearance);
    }

    TEST(Clearance_Basic_P1_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.003, clearance);
    }

    TEST(Clearance_Reinforced_P1_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0055, clearance);
    }

    TEST(Clearance_Basic_P2_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.003, clearance);
    }

    TEST(Clearance_Reinforced_P2_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0055, clearance);
    }

    TEST(Clearance_Basic_P3_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.003, clearance);
    }

    TEST(Clearance_Reinforced_P3_OVC_III_Low_Altitude_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_III;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0055, clearance);
    }

    TEST(Clearance_Basic_P1_OVC_I_High_Altitude_Low_Frequency) {
        altitude.set_maximum(5000);
        mainSupplyVoltage.set_nominal(250);
        frequency = 30000;
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.00074, clearance);
    }

    TEST(Clearance_Printed_Basic) {
        maximumVoltageRms = 120;
        maximumVoltagePeak = 200;
        altitude.set_maximum(2000);
        frequency = 30000;
        maximumVoltagePeak = 2000;
        mainSupplyVoltage.set_nominal(120);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::PRINTED);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.0002, clearance);
    }

    TEST(Clearance_Basic_P1_OVC_I_Low_Altitude_High_Frequency) {
        mainSupplyVoltage.set_nominal(250);
        frequency = 500000;
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto clearance = standard.calculate_clearance(inputs);
        CHECK_EQUAL(0.011, clearance);
    }

}


SUITE(CreepageDistance_IEC_60335) {

    auto standard = OpenMagnetics::InsulationIEC60335Model();
    auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_603351};
    auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_I;
    double maximumVoltageRms = 250;
    double maximumVoltagePeak = 400;
    double frequency = 30000;
    OpenMagnetics::DimensionWithTolerance altitude;
    OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;


    TEST(Creepage_Distance_Functional_P1_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::FUNCTIONAL;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.00075, creepageDistance);
    }
    TEST(Creepage_Distance_Basic_P1_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.001, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.002, creepageDistance);
    }

    TEST(Creepage_Distance_Functional_P2_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::FUNCTIONAL;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0016, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.002, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.004, creepageDistance);
    }

    TEST(Creepage_Distance_Functional_P3_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::FUNCTIONAL;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.004, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.005, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_I_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.01, creepageDistance);
    }

    TEST(Creepage_Distance_Functional_P1_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::FUNCTIONAL;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.00075, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.001, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.002, creepageDistance);
    }

    TEST(Creepage_Distance_Functional_P2_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::FUNCTIONAL;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0022, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0028, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0056, creepageDistance);
    }

    TEST(Creepage_Distance_Functional_P3_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::FUNCTIONAL;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0045, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0056, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_II_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0112, creepageDistance);
    }

    TEST(Creepage_Distance_Functional_P1_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::FUNCTIONAL;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.00075, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.001, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.002, creepageDistance);
    }

    TEST(Creepage_Distance_Functional_P2_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::FUNCTIONAL;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0032, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P2_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.004, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P2_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.008, creepageDistance);
    }

    TEST(Creepage_Distance_Functional_P3_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::FUNCTIONAL;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.005, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P3_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0063, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P3_GROUP_IIIA_Low_Frequency) {
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_IIIA;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P3;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.01260, creepageDistance);
    }

    TEST(Creepage_Distance_Functional_P1_GROUP_I_High_Frequency) {
        frequency = 700000;
        maximumVoltagePeak = 800;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::FUNCTIONAL;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0038, creepageDistance);
    }

    TEST(Creepage_Distance_Basic_P1_GROUP_I_High_Frequency) {
        frequency = 700000;
        maximumVoltagePeak = 800;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0038, creepageDistance);
    }

    TEST(Creepage_Distance_Reinforced_P1_GROUP_I_High_Frequency) {
        frequency = 700000;
        maximumVoltagePeak = 800;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(250);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK_EQUAL(0.0076, creepageDistance);
    }
}