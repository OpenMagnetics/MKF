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
        auto leftWire = OpenMagnetics::find_wire_by_name("Litz Round 0.016 - Grade 1");
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