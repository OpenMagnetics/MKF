#include "Insulation.h"
#include "Painter.h"
#include "CoreAdviser.h"
#include "CoilAdviser.h"
#include "InputsWrapper.h"
#include "TestingUtils.h"

#include <UnitTest++.h>
#include <vector>

SUITE(SolidInsulationRequirements) {

    auto standardCoordinator = OpenMagnetics::InsulationCoordinator();
    auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
    auto cti = OpenMagnetics::Cti::GROUP_I;
    double maximumVoltageRms = 666;
    double maximumVoltagePeak = 800;
    double frequency = 30000;
    OpenMagnetics::DimensionWithTolerance altitude;
    OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;
    auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;


    TEST(Test_CoilAdviser_Wires_Withstand_Voltage_Same_Isolation_Sides_Basic_No_FIW) {

        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(1);
        inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
        inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<OpenMagnetics::IsolationSide>{OpenMagnetics::IsolationSide::PRIMARY, OpenMagnetics::IsolationSide::PRIMARY});

        OpenMagnetics::CoilAdviser coilAdviser;
        auto withstandVoltageForWires = coilAdviser.get_solid_insulation_requirements_for_wires(inputs);
        CHECK_EQUAL(2, withstandVoltageForWires.size());
        CHECK_EQUAL(2, withstandVoltageForWires[0].size());
        CHECK(withstandVoltageForWires[0][0].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[0][1].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[0][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[0][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[0][0].get_minimum_grade());

        CHECK(withstandVoltageForWires[0][0].get_minimum_grade().value() == 1);
        CHECK(withstandVoltageForWires[1][0].get_minimum_breakdown_voltage() > 0);
        CHECK(withstandVoltageForWires[1][1].get_minimum_breakdown_voltage() > 0);
        CHECK(withstandVoltageForWires[1][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[1][0].get_minimum_number_layers().value() == 1);
        CHECK(!withstandVoltageForWires[1][0].get_minimum_grade());
    }

    TEST(Test_CoilAdviser_Wires_Withstand_Voltage_Different_Isolation_Sides_Functional_No_FIW) {
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto insulationType = OpenMagnetics::InsulationType::FUNCTIONAL;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(1);
        inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
        inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<OpenMagnetics::IsolationSide>{OpenMagnetics::IsolationSide::PRIMARY, OpenMagnetics::IsolationSide::SECONDARY});

        OpenMagnetics::CoilAdviser coilAdviser;
        auto withstandVoltageForWires = coilAdviser.get_solid_insulation_requirements_for_wires(inputs);
        CHECK_EQUAL(1, withstandVoltageForWires.size());
        CHECK_EQUAL(2, withstandVoltageForWires[0].size());

        CHECK(withstandVoltageForWires[0][0].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[0][1].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[0][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[0][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[0][0].get_minimum_grade());
        CHECK(withstandVoltageForWires[0][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[0][1].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[0][1].get_minimum_grade());
    }

    TEST(Test_CoilAdviser_Wires_Withstand_Voltage_Different_Isolation_Sides_Basic_No_FIW) {
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(1);
        inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
        inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<OpenMagnetics::IsolationSide>{OpenMagnetics::IsolationSide::PRIMARY, OpenMagnetics::IsolationSide::SECONDARY});

        OpenMagnetics::CoilAdviser coilAdviser;
        auto withstandVoltageForWires = coilAdviser.get_solid_insulation_requirements_for_wires(inputs);
        CHECK_EQUAL(3, withstandVoltageForWires.size());
        CHECK_EQUAL(2, withstandVoltageForWires[0].size());

        CHECK(withstandVoltageForWires[0][0].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[0][1].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[0][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[0][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[0][0].get_minimum_grade());
        CHECK(withstandVoltageForWires[0][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[0][1].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[0][1].get_minimum_grade());

        CHECK(withstandVoltageForWires[1][0].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[1][1].get_minimum_breakdown_voltage() == 6000);
        CHECK(withstandVoltageForWires[1][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[1][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[1][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[1][1].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[1][0].get_minimum_grade());
        CHECK(!withstandVoltageForWires[1][1].get_minimum_grade());

        CHECK(withstandVoltageForWires[2][0].get_minimum_breakdown_voltage() == 6000);
        CHECK(withstandVoltageForWires[2][1].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[2][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[2][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[2][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[2][1].get_minimum_number_layers().value() == 1);
        CHECK(!withstandVoltageForWires[2][0].get_minimum_grade());
        CHECK(withstandVoltageForWires[2][1].get_minimum_grade());
    }

    TEST(Test_CoilAdviser_Wires_Withstand_Voltage_Different_Isolation_Sides_Supplementary_No_FIW) {
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto insulationType = OpenMagnetics::InsulationType::SUPPLEMENTARY;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(1);
        inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
        inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<OpenMagnetics::IsolationSide>{OpenMagnetics::IsolationSide::PRIMARY, OpenMagnetics::IsolationSide::SECONDARY});

        OpenMagnetics::CoilAdviser coilAdviser;
        auto withstandVoltageForWires = coilAdviser.get_solid_insulation_requirements_for_wires(inputs);
        CHECK_EQUAL(3, withstandVoltageForWires.size());
        CHECK_EQUAL(2, withstandVoltageForWires[0].size());

        CHECK(withstandVoltageForWires[0][0].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[0][1].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[0][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[0][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[0][0].get_minimum_grade());
        CHECK(withstandVoltageForWires[0][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[0][1].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[0][1].get_minimum_grade());

        CHECK(withstandVoltageForWires[1][0].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[1][1].get_minimum_breakdown_voltage() == 6000);
        CHECK(withstandVoltageForWires[1][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[1][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[1][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[1][1].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[1][0].get_minimum_grade());
        CHECK(!withstandVoltageForWires[1][1].get_minimum_grade());

        CHECK(withstandVoltageForWires[2][0].get_minimum_breakdown_voltage() == 6000);
        CHECK(withstandVoltageForWires[2][1].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[2][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[2][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[2][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[2][1].get_minimum_number_layers().value() == 1);
        CHECK(!withstandVoltageForWires[2][0].get_minimum_grade());
        CHECK(withstandVoltageForWires[2][1].get_minimum_grade());
    }

    TEST(Test_CoilAdviser_Wires_Withstand_Voltage_Different_Isolation_Sides_Reinforced_No_FIW) {
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(1);
        inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
        inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<OpenMagnetics::IsolationSide>{OpenMagnetics::IsolationSide::PRIMARY, OpenMagnetics::IsolationSide::SECONDARY});

        OpenMagnetics::CoilAdviser coilAdviser;
        auto withstandVoltageForWires = coilAdviser.get_solid_insulation_requirements_for_wires(inputs);
        CHECK_EQUAL(3, withstandVoltageForWires.size());
        CHECK_EQUAL(2, withstandVoltageForWires[0].size());

        CHECK(withstandVoltageForWires[0][0].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[0][1].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[0][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[0][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[0][0].get_minimum_grade());
        CHECK(withstandVoltageForWires[0][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[0][1].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[0][1].get_minimum_grade());

        CHECK(withstandVoltageForWires[1][0].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[1][1].get_minimum_breakdown_voltage() == 8000);
        CHECK(withstandVoltageForWires[1][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[1][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[1][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[1][1].get_minimum_number_layers().value() == 3);
        CHECK(withstandVoltageForWires[1][0].get_minimum_grade());
        CHECK(!withstandVoltageForWires[1][1].get_minimum_grade());

        CHECK(withstandVoltageForWires[2][0].get_minimum_breakdown_voltage() == 8000);
        CHECK(withstandVoltageForWires[2][1].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[2][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[2][0].get_minimum_number_layers().value() == 3);
        CHECK(withstandVoltageForWires[2][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[2][1].get_minimum_number_layers().value() == 1);
        CHECK(!withstandVoltageForWires[2][0].get_minimum_grade());
        CHECK(withstandVoltageForWires[2][1].get_minimum_grade());
    }

    TEST(Test_CoilAdviser_Wires_Withstand_Voltage_Same_Isolation_Sides_Basic_FIW) {

        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_623681};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(1);
        inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
        inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<OpenMagnetics::IsolationSide>{OpenMagnetics::IsolationSide::PRIMARY, OpenMagnetics::IsolationSide::PRIMARY});

        OpenMagnetics::CoilAdviser coilAdviser;
        auto withstandVoltageForWires = coilAdviser.get_solid_insulation_requirements_for_wires(inputs);
        CHECK_EQUAL(2, withstandVoltageForWires.size());
        CHECK_EQUAL(2, withstandVoltageForWires[0].size());
        CHECK(withstandVoltageForWires[0][0].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[0][1].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[0][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[0][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[0][0].get_minimum_grade());

        CHECK(withstandVoltageForWires[0][0].get_minimum_grade().value() == 1);
        CHECK(withstandVoltageForWires[1][0].get_minimum_breakdown_voltage() > 0);
        CHECK(withstandVoltageForWires[1][1].get_minimum_breakdown_voltage() > 0);
        CHECK(withstandVoltageForWires[1][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[1][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[1][0].get_minimum_grade());
        CHECK(withstandVoltageForWires[1][0].get_minimum_grade().value() == 3);
    }

    TEST(Test_CoilAdviser_Wires_Withstand_Voltage_Different_Isolation_Sides_Functional_FIW) {
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_623681};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto insulationType = OpenMagnetics::InsulationType::FUNCTIONAL;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(1);
        inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
        inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<OpenMagnetics::IsolationSide>{OpenMagnetics::IsolationSide::PRIMARY, OpenMagnetics::IsolationSide::SECONDARY});

        OpenMagnetics::CoilAdviser coilAdviser;
        auto withstandVoltageForWires = coilAdviser.get_solid_insulation_requirements_for_wires(inputs);
        CHECK_EQUAL(1, withstandVoltageForWires.size());
        CHECK_EQUAL(2, withstandVoltageForWires[0].size());

        CHECK(withstandVoltageForWires[0][0].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[0][1].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[0][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[0][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[0][0].get_minimum_grade());
        CHECK(withstandVoltageForWires[0][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[0][1].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[0][1].get_minimum_grade());
    }

    TEST(Test_CoilAdviser_Wires_Withstand_Voltage_Different_Isolation_Sides_Basic_FIW) {
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_623681};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(1);
        inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
        inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<OpenMagnetics::IsolationSide>{OpenMagnetics::IsolationSide::PRIMARY, OpenMagnetics::IsolationSide::SECONDARY});

        OpenMagnetics::CoilAdviser coilAdviser;
        auto withstandVoltageForWires = coilAdviser.get_solid_insulation_requirements_for_wires(inputs);
        CHECK_EQUAL(3, withstandVoltageForWires.size());
        CHECK_EQUAL(2, withstandVoltageForWires[0].size());

        CHECK(withstandVoltageForWires[0][0].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[0][1].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[0][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[0][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[0][0].get_minimum_grade());
        CHECK(withstandVoltageForWires[0][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[0][1].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[0][1].get_minimum_grade());

        CHECK(withstandVoltageForWires[1][0].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[1][1].get_minimum_breakdown_voltage() == 2500);
        CHECK(withstandVoltageForWires[1][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[1][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[1][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[1][1].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[1][0].get_minimum_grade());
        CHECK(withstandVoltageForWires[1][1].get_minimum_grade());
        CHECK(withstandVoltageForWires[1][1].get_minimum_grade().value() == 3);

        CHECK(withstandVoltageForWires[2][0].get_minimum_breakdown_voltage() == 2500);
        CHECK(withstandVoltageForWires[2][1].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[2][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[2][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[2][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[2][1].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[2][0].get_minimum_grade());
        CHECK(withstandVoltageForWires[2][0].get_minimum_grade().value() == 3);
        CHECK(withstandVoltageForWires[2][1].get_minimum_grade());
    }

    TEST(Test_CoilAdviser_Wires_Withstand_Voltage_Different_Isolation_Sides_Supplementary_FIW) {
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_623681};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto insulationType = OpenMagnetics::InsulationType::SUPPLEMENTARY;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(1);
        inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
        inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<OpenMagnetics::IsolationSide>{OpenMagnetics::IsolationSide::PRIMARY, OpenMagnetics::IsolationSide::SECONDARY});

        OpenMagnetics::CoilAdviser coilAdviser;
        auto withstandVoltageForWires = coilAdviser.get_solid_insulation_requirements_for_wires(inputs);
        CHECK_EQUAL(3, withstandVoltageForWires.size());
        CHECK_EQUAL(2, withstandVoltageForWires[0].size());

        CHECK(withstandVoltageForWires[0][0].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[0][1].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[0][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[0][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[0][0].get_minimum_grade());
        CHECK(withstandVoltageForWires[0][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[0][1].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[0][1].get_minimum_grade());

        CHECK(withstandVoltageForWires[1][0].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[1][1].get_minimum_breakdown_voltage() == 2500);
        CHECK(withstandVoltageForWires[1][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[1][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[1][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[1][1].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[1][0].get_minimum_grade());
        CHECK(withstandVoltageForWires[1][1].get_minimum_grade());
        CHECK(withstandVoltageForWires[1][1].get_minimum_grade().value() == 3);

        CHECK(withstandVoltageForWires[2][0].get_minimum_breakdown_voltage() == 2500);
        CHECK(withstandVoltageForWires[2][1].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[2][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[2][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[2][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[2][1].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[2][0].get_minimum_grade());
        CHECK(withstandVoltageForWires[2][0].get_minimum_grade().value() == 3);
        CHECK(withstandVoltageForWires[2][1].get_minimum_grade());
    }

    TEST(Test_CoilAdviser_Wires_Withstand_Voltage_Different_Isolation_Sides_Reinforced_FIW) {
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_623681};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(1);
        inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
        inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<OpenMagnetics::IsolationSide>{OpenMagnetics::IsolationSide::PRIMARY, OpenMagnetics::IsolationSide::SECONDARY});

        OpenMagnetics::CoilAdviser coilAdviser;
        auto withstandVoltageForWires = coilAdviser.get_solid_insulation_requirements_for_wires(inputs);
        CHECK_EQUAL(3, withstandVoltageForWires.size());
        CHECK_EQUAL(2, withstandVoltageForWires[0].size());

        CHECK(withstandVoltageForWires[0][0].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[0][1].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[0][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[0][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[0][0].get_minimum_grade());
        CHECK(withstandVoltageForWires[0][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[0][1].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[0][1].get_minimum_grade()); 

        CHECK(withstandVoltageForWires[1][0].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[1][1].get_minimum_breakdown_voltage() == 5000);
        CHECK(withstandVoltageForWires[1][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[1][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[1][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[1][1].get_minimum_number_layers().value() == 3);
        CHECK(withstandVoltageForWires[1][0].get_minimum_grade());
        CHECK(withstandVoltageForWires[1][1].get_minimum_grade());
        CHECK(withstandVoltageForWires[1][1].get_minimum_grade().value() == 3);

        CHECK(withstandVoltageForWires[2][0].get_minimum_breakdown_voltage() == 5000);
        CHECK(withstandVoltageForWires[2][1].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[2][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[2][0].get_minimum_number_layers().value() == 3);
        CHECK(withstandVoltageForWires[2][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[2][1].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[2][0].get_minimum_grade());
        CHECK(withstandVoltageForWires[2][0].get_minimum_grade().value() == 3);
        CHECK(withstandVoltageForWires[2][1].get_minimum_grade());
    }

    TEST(Test_CoilAdviser_Wires_Withstand_Voltage_Different_Isolation_Sides_Double_FIW) {
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_623681};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto insulationType = OpenMagnetics::InsulationType::DOUBLE;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(1);
        inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
        inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<OpenMagnetics::IsolationSide>{OpenMagnetics::IsolationSide::PRIMARY, OpenMagnetics::IsolationSide::SECONDARY});

        OpenMagnetics::CoilAdviser coilAdviser;
        auto withstandVoltageForWires = coilAdviser.get_solid_insulation_requirements_for_wires(inputs);
        CHECK_EQUAL(6, withstandVoltageForWires.size());
        CHECK_EQUAL(2, withstandVoltageForWires[0].size());

        CHECK(withstandVoltageForWires[0][0].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[0][1].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[0][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[0][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[0][0].get_minimum_grade());
        CHECK(withstandVoltageForWires[0][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[0][1].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[0][1].get_minimum_grade());

        CHECK(withstandVoltageForWires[1][0].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[1][1].get_minimum_breakdown_voltage() == 5000);
        CHECK(withstandVoltageForWires[1][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[1][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[1][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[1][1].get_minimum_number_layers().value() == 3);
        CHECK(withstandVoltageForWires[1][0].get_minimum_grade());
        CHECK(withstandVoltageForWires[1][1].get_minimum_grade());
        CHECK(withstandVoltageForWires[1][1].get_minimum_grade().value() == 3);

        CHECK(withstandVoltageForWires[2][0].get_minimum_breakdown_voltage() == 5000);
        CHECK(withstandVoltageForWires[2][1].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[2][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[2][0].get_minimum_number_layers().value() == 3);
        CHECK(withstandVoltageForWires[2][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[2][1].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[2][0].get_minimum_grade());
        CHECK(withstandVoltageForWires[2][0].get_minimum_grade().value() == 3);
        CHECK(withstandVoltageForWires[2][1].get_minimum_grade());

        CHECK(withstandVoltageForWires[3][0].get_minimum_breakdown_voltage() == 2500);
        CHECK(withstandVoltageForWires[3][1].get_minimum_breakdown_voltage() == 2500);
        CHECK(withstandVoltageForWires[3][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[3][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[3][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[3][1].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[3][0].get_minimum_grade());
        CHECK(withstandVoltageForWires[3][0].get_minimum_grade().value() == 3);
        CHECK(withstandVoltageForWires[3][1].get_minimum_grade());
        CHECK(withstandVoltageForWires[3][1].get_minimum_grade().value() == 3);

        CHECK(withstandVoltageForWires[4][0].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[4][1].get_minimum_breakdown_voltage() == 2500);
        CHECK(withstandVoltageForWires[4][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[4][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[4][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[4][1].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[4][0].get_minimum_grade());
        CHECK(withstandVoltageForWires[4][1].get_minimum_grade());
        CHECK(withstandVoltageForWires[4][1].get_minimum_grade().value() == 3);

        CHECK(withstandVoltageForWires[5][0].get_minimum_breakdown_voltage() == 2500);
        CHECK(withstandVoltageForWires[5][1].get_minimum_breakdown_voltage() == 0);
        CHECK(withstandVoltageForWires[5][0].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[5][0].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[5][1].get_minimum_number_layers());
        CHECK(withstandVoltageForWires[5][1].get_minimum_number_layers().value() == 1);
        CHECK(withstandVoltageForWires[5][0].get_minimum_grade());
        CHECK(withstandVoltageForWires[5][0].get_minimum_grade().value() == 3);
        CHECK(withstandVoltageForWires[5][1].get_minimum_grade());
    }

    // TEST(Test_CoilAdviser_Wires_Withstand_Voltage_Different_Isolation_Sides_Reinforced) {
    //     auto standardCoordinator = OpenMagnetics::InsulationCoordinator();
    //     auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
    //     auto cti = OpenMagnetics::Cti::GROUP_I;
    //     double maximumVoltageRms = 666;
    //     double maximumVoltagePeak = 800;
    //     double frequency = 30000;
    //     OpenMagnetics::DimensionWithTolerance altitude;
    //     OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;

    //     auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641};
    //     altitude.set_maximum(2000);
    //     mainSupplyVoltage.set_nominal(400);
    //     auto insulationType = OpenMagnetics::InsulationType::REINFORCED;
    //     auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
    //     OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
    //     OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
    //     dimensionWithTolerance.set_nominal(1);
    //     inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
    //     inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<OpenMagnetics::IsolationSide>{OpenMagnetics::IsolationSide::PRIMARY, OpenMagnetics::IsolationSide::SECONDARY});

    //     OpenMagnetics::CoilAdviser coilAdviser;
    //     auto withstandVoltageForWires = coilAdviser.get_withstand_voltage_for_wires(inputs);
    //     CHECK_EQUAL(1, withstandVoltageForWires.size());
    //     CHECK_EQUAL(2, withstandVoltageForWires[0].size());
    //     CHECK(withstandVoltageForWires[0][0] > 0);
    //     CHECK(withstandVoltageForWires[0][1] > 0);
    //     std::cout << withstandVoltageForWires[0][1] << std::endl;
    // }

    // TEST(Test_CoilAdviser_Wires_Withstand_Voltage_Three_Different_Isolation_Sides) {
    //     auto standardCoordinator = OpenMagnetics::InsulationCoordinator();
    //     auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
    //     auto cti = OpenMagnetics::Cti::GROUP_I;
    //     double maximumVoltageRms = 666;
    //     double maximumVoltagePeak = 800;
    //     double frequency = 30000;
    //     OpenMagnetics::DimensionWithTolerance altitude;
    //     OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;

    //     auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641};
    //     altitude.set_maximum(2000);
    //     mainSupplyVoltage.set_nominal(400);
    //     auto insulationType = OpenMagnetics::InsulationType::BASIC;
    //     auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
    //     OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
    //     OpenMagnetics::DimensionWithTolerance dimensionWithTolerance;
    //     dimensionWithTolerance.set_nominal(1);
    //     inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance, dimensionWithTolerance});
    //     inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<OpenMagnetics::IsolationSide>{OpenMagnetics::IsolationSide::PRIMARY, OpenMagnetics::IsolationSide::SECONDARY, OpenMagnetics::IsolationSide::TERTIARY});

    //     OpenMagnetics::CoilAdviser coilAdviser;
    //     auto withstandVoltageForWires = coilAdviser.get_withstand_voltage_for_wires(inputs);
    //     CHECK_EQUAL(3, withstandVoltageForWires.size());
    //     CHECK_EQUAL(3, withstandVoltageForWires[0].size());
    //     CHECK(withstandVoltageForWires[0][0] == 0);
    //     CHECK(withstandVoltageForWires[0][1] > 0);
    //     CHECK(withstandVoltageForWires[0][2] > 0);
    //     CHECK(withstandVoltageForWires[1][0] > 0);
    //     CHECK(withstandVoltageForWires[1][1] == 0);
    //     CHECK(withstandVoltageForWires[1][2] > 0);
    //     CHECK(withstandVoltageForWires[2][0] > 0);
    //     CHECK(withstandVoltageForWires[2][1] > 0);
    //     CHECK(withstandVoltageForWires[2][2] == 0);
    // }
}

SUITE(CoilAdviser) {
    TEST(Test_CoilAdviser_Random) {
        srand (time(NULL));

        int count = 10;
        while (count > 0) {
            std::vector<double> turnsRatios;
            int64_t numberStacks = 1;

            std::vector<int64_t> numberTurns;
            std::vector<int64_t> numberParallels;
            std::vector<OpenMagnetics::IsolationSide> isolationSides;
            int64_t numberPhysicalTurns = std::numeric_limits<int64_t>::max();
            // for (size_t windingIndex = 0; windingIndex < std::rand() % 2 + 1UL; ++windingIndex)
            for (size_t windingIndex = 0; windingIndex < std::rand() % 4 + 1UL; ++windingIndex)
            {
                numberTurns.push_back(std::rand() % 100 + 1L);
                isolationSides.push_back(OpenMagnetics::get_isolation_side_from_index(static_cast<size_t>(std::rand() % 10 + 1L)));
                numberParallels.push_back(std::rand() % 3 + 1L);
                numberPhysicalTurns = std::min(numberPhysicalTurns, numberTurns.back() * numberParallels.back());
            }
            for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
                turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
            }

            double frequency = std::rand() % 1000000 + 10000L;
            double magnetizingInductance = double(std::rand() % 10000) *  1e-6;
            double temperature = 25;
            OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::SINUSOIDAL;
            double peakToPeak = std::rand() % 30 + 1L;
            double dutyCycle = double(std::rand() % 100) / 100;
            double dcCurrent = 0;
            if (numberTurns.size() > 1) {
                dcCurrent = std::rand() % 30;
            }
            auto coreShapeNames = OpenMagnetics::get_shape_names(false);
            std::string coreShapeName;
            OpenMagnetics::MagneticWrapper magnetic;

            double gapLength = 1;
            double columnHeight = 1;
            int numberGaps = 1;
            int gappingTypeIndex = std::rand() % 4;
            OpenMagnetics::GappingType gappingType = magic_enum::enum_cast<OpenMagnetics::GappingType>(gappingTypeIndex).value();
            if (gappingType == OpenMagnetics::GappingType::DISTRIBUTED) {
                numberGaps = std::rand() % 3;
                numberGaps *= 2;
                numberGaps += 3;
            }

            while (true) {
                coreShapeName = coreShapeNames[rand() % coreShapeNames.size()];
                auto shape = OpenMagnetics::find_core_shape_by_name(coreShapeName);
                if (shape.get_family() == OpenMagnetics::CoreShapeFamily::PQI || 
                    shape.get_family() == OpenMagnetics::CoreShapeFamily::UI || 
                    shape.get_family() == OpenMagnetics::CoreShapeFamily::EI) {
                    continue;
                }

                do {
                    auto shape = OpenMagnetics::find_core_shape_by_name(coreShapeName);
                    auto corePiece = OpenMagnetics::CorePiece::factory(shape);
                    gapLength = double(rand() % 10000 + 1.0) / 1000000;
                    columnHeight = corePiece->get_columns()[0].get_height();
                }
                while (columnHeight < (gapLength * numberGaps));

                std::vector<OpenMagnetics::CoreGap> gapping;
                switch (gappingType) {
                    case OpenMagnetics::GappingType::GRINDED:
                        gapping = OpenMagneticsTesting::get_grinded_gap(gapLength);
                        break;
                    case OpenMagnetics::GappingType::SPACER:
                        gapping = OpenMagneticsTesting::get_spacer_gap(gapLength);
                        break;
                    case OpenMagnetics::GappingType::RESIDUAL:
                        gapping = OpenMagneticsTesting::get_residual_gap();
                        break;
                    case OpenMagnetics::GappingType::DISTRIBUTED:
                        gapping = OpenMagneticsTesting::get_distributed_gap(gapLength, numberGaps);
                        break;
                }
                
                magnetic = OpenMagneticsTesting::get_quick_magnetic(coreShapeName, gapping, numberTurns, numberStacks, "3C91");
                break;
            }

            int turnsAlignmentIndex = std::rand() % 4;
            int windingOrientationIndex = (std::rand() % 2) * 2;  // To avoid 1, which is RADIAL
            int layersOrientationIndex = (std::rand() % 2) * 2;  // To avoid 1, which is RADIAL
            OpenMagnetics::CoilAlignment turnsAlignment = magic_enum::enum_cast<OpenMagnetics::CoilAlignment>(turnsAlignmentIndex).value();
            OpenMagnetics::WindingOrientation windingOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(windingOrientationIndex).value();
            OpenMagnetics::WindingOrientation layersOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(layersOrientationIndex).value();


            for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
                std::cout << "numberTurns: " << numberTurns[windingIndex] << std::endl;
            }
            for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
                std::cout << "isolationSides: " << magic_enum::enum_name(isolationSides[windingIndex]) << std::endl;
            }
            std::cout << "frequency: " << frequency << std::endl;
            std::cout << "peakToPeak: " << peakToPeak << std::endl;
            std::cout << "magnetizingInductance: " << magnetizingInductance << std::endl;
            std::cout << "dutyCycle: " << dutyCycle << std::endl;
            std::cout << "dcCurrent: " << dcCurrent << std::endl;
            std::cout << "coreShapeName: " << coreShapeName << std::endl;
            std::cout << "gappingTypeIndex: " << gappingTypeIndex << std::endl;
            std::cout << "gapLength: " << gapLength << std::endl;
            std::cout << "numberGaps: " << numberGaps << std::endl;
            std::cout << "magnetic.get_mutable_core().get_shape_family(): " << magic_enum::enum_name(magnetic.get_mutable_core().get_shape_family()) << std::endl;
            std::cout << "turnsAlignmentIndex: " << turnsAlignmentIndex << std::endl;
            std::cout << "windingOrientationIndex: " << windingOrientationIndex << std::endl;
            std::cout << "layersOrientationIndex: " << layersOrientationIndex << std::endl;

            magnetic.get_mutable_coil().set_turns_alignment(turnsAlignment);
            magnetic.get_mutable_coil().set_winding_orientation(windingOrientation);
            magnetic.get_mutable_coil().set_layers_orientation(layersOrientation);

            auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  waveShape,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  dcCurrent,
                                                                                                  turnsRatios);

            inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
            OpenMagnetics::MasWrapper masMagnetic;
            inputs.process_waveforms();
            masMagnetic.set_inputs(inputs);
            masMagnetic.set_magnetic(magnetic);

            OpenMagnetics::CoilAdviser coilAdviser;
            auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

            if (masMagneticsWithCoil.size() > 0) {
                auto masMagneticWithCoil = masMagneticsWithCoil[0].first;
                if (!masMagneticWithCoil.get_magnetic().get_coil().get_turns_description()) {
                    continue;
                }
                OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
                std::string filePath = __FILE__;
                auto randNumber =  std::rand();
                // auto outputFilePath = filePath.substr(0, filePath.rfind("/")).append("/../output/");
                // {
                //     auto outFile = outputFilePath;
                //     std::string filename = "Test_CoilAdviser" + std::to_string(randNumber) + "_Quiver.svg";
                //     outFile.append(filename);
                //     OpenMagnetics::Painter painter(outFile, OpenMagnetics::Painter::PainterModes::QUIVER);
                //     painter.set_number_points_x(50);
                //     painter.set_number_points_y(50);
                //     painter.set_logarithmic_scale(false);
                //     painter.set_fringing_effect(true);
                //     painter.set_maximum_scale_value(std::nullopt);
                //     painter.set_minimum_scale_value(std::nullopt);
                //     painter.paint_magnetic_field(inputs.get_operating_point(0), masMagneticWithCoil.get_magnetic());
                //     painter.paint_core(masMagneticWithCoil.get_magnetic());
                //     painter.paint_bobbin(masMagneticWithCoil.get_magnetic());
                //     painter.paint_coil_turns(masMagneticWithCoil.get_magnetic());
                //     painter.export_svg();
                // }
                // {
                //     auto outFile = outputFilePath;
                //     std::string filename = "Test_CoilAdviser" + std::to_string(randNumber) + "_Quiver_Log.svg";
                //     outFile.append(filename);
                //     OpenMagnetics::Painter painter(outFile, OpenMagnetics::Painter::PainterModes::QUIVER);
                //     painter.set_number_points_x(50);
                //     painter.set_number_points_y(50);
                //     painter.set_logarithmic_scale(true);
                //     painter.set_fringing_effect(true);
                //     painter.set_maximum_scale_value(std::nullopt);
                //     painter.set_minimum_scale_value(std::nullopt);
                //     painter.paint_magnetic_field(inputs.get_operating_point(0), masMagneticWithCoil.get_magnetic());
                //     painter.paint_core(masMagneticWithCoil.get_magnetic());
                //     painter.paint_bobbin(masMagneticWithCoil.get_magnetic());
                //     painter.paint_coil_turns(masMagneticWithCoil.get_magnetic());
                //     painter.export_svg();
                // }
                // {
                //     auto outFile = outputFilePath;
                //     std::string filename = "Test_CoilAdviser" + std::to_string(randNumber) + "_Contour.svg";
                //     outFile.append(filename);
                //     OpenMagnetics::Painter painter(outFile, OpenMagnetics::Painter::PainterModes::CONTOUR);
                //     painter.set_number_points_x(50);
                //     painter.set_number_points_y(50);
                //     painter.set_logarithmic_scale(false);
                //     painter.set_fringing_effect(true);
                //     painter.set_maximum_scale_value(std::nullopt);
                //     painter.set_minimum_scale_value(std::nullopt);
                //     painter.paint_magnetic_field(inputs.get_operating_point(0), masMagneticWithCoil.get_magnetic());
                //     painter.paint_core(masMagneticWithCoil.get_magnetic());
                //     painter.paint_bobbin(masMagneticWithCoil.get_magnetic());
                //     painter.paint_coil_turns(masMagneticWithCoil.get_magnetic());
                //     painter.export_svg();
                // }
                // {
                //     auto outFile = outputFilePath;
                //     std::string filename = "Test_CoilAdviser" + std::to_string(randNumber) + "_Contour_Log.svg";
                //     outFile.append(filename);
                //     OpenMagnetics::Painter painter(outFile, OpenMagnetics::Painter::PainterModes::CONTOUR);
                //     painter.set_number_points_x(50);
                //     painter.set_number_points_y(50);
                //     painter.set_logarithmic_scale(true);
                //     painter.set_fringing_effect(true);
                //     painter.set_maximum_scale_value(std::nullopt);
                //     painter.set_minimum_scale_value(std::nullopt);
                //     painter.paint_magnetic_field(inputs.get_operating_point(0), masMagneticWithCoil.get_magnetic());
                //     painter.paint_core(masMagneticWithCoil.get_magnetic());
                //     painter.paint_bobbin(masMagneticWithCoil.get_magnetic());
                //     painter.paint_coil_turns(masMagneticWithCoil.get_magnetic());
                //     painter.export_svg();
                // }
                count--;
            }
        }

    }

    TEST(Test_CoilAdviser_Random_0) {
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.003);
        std::vector<double> turnsRatios;
        int64_t numberStacks = 1;


        double magnetizingInductance = 10e-6;
        double temperature = 25;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double dutyCycle = 0.5;
        double dcCurrent = 0;

        std::vector<int64_t> numberTurns = {82, 5};
        std::vector<int64_t> numberParallels = {1, 1};
        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }
        double frequency = 675590;
        double peakToPeak = 26;

        auto magnetic = OpenMagneticsTesting::get_quick_magnetic("EP 20", gapping, numberTurns, numberStacks, "3C91");
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        OpenMagnetics::MasWrapper masMagnetic;
        inputs.process_waveforms();
        masMagnetic.set_inputs(inputs);
        masMagnetic.set_magnetic(magnetic);

        OpenMagnetics::CoilAdviser coilAdviser;
        auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

        if (masMagneticsWithCoil.size() > 0) {
            auto masMagneticWithCoil = masMagneticsWithCoil[0].first;
            OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
        }
    }
    TEST(Test_CoilAdviser_Random_1) {
        std::vector<double> turnsRatios;
        int64_t numberStacks = 1;


        double magnetizingInductance = 10e-6;
        double temperature = 25;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double dutyCycle = 0.5;
        double dcCurrent = 0;

        std::vector<int64_t> numberTurns = {16, 34};
        std::vector<int64_t> numberParallels = {1, 1};
        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }
        double frequency = 811022;
        double peakToPeak = 1;

        OpenMagnetics::MagneticWrapper magnetic;
        std::string coreShapeName  = "P 7.35X3.6";
        double gapLength = 1;
        double columnHeight = 1;
        do {
            auto shape = OpenMagnetics::find_core_shape_by_name(coreShapeName);
            auto corePiece = OpenMagnetics::CorePiece::factory(shape);
            gapLength = double(rand() % 10000 + 1.0) / 1000000;
            columnHeight = corePiece->get_columns()[0].get_height();
        }
        while (columnHeight < gapLength);
        auto gapping = OpenMagneticsTesting::get_grinded_gap(gapLength);
        magnetic = OpenMagneticsTesting::get_quick_magnetic(coreShapeName, gapping, numberTurns, numberStacks, "3C91");
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        OpenMagnetics::MasWrapper masMagnetic;
        inputs.process_waveforms();
        masMagnetic.set_inputs(inputs);
        masMagnetic.set_magnetic(magnetic);

        OpenMagnetics::CoilAdviser coilAdviser;
        auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

        if (masMagneticsWithCoil.size() > 0) {
            auto masMagneticWithCoil = masMagneticsWithCoil[0].first;
            OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
            std::string filePath = __FILE__;
            auto outputFilePath = filePath.substr(0, filePath.rfind("/")).append("/../output/");
            auto outFile = outputFilePath;
            std::string filename = "Test_CoilAdviser" + std::to_string(std::rand()) + ".svg";
            outFile.append(filename);
            // OpenMagnetics::Painter painter(outFile);

            // painter.paint_core(masMagneticWithCoil.get_magnetic());
            // painter.paint_bobbin(masMagneticWithCoil.get_magnetic());
            // painter.paint_coil_turns(masMagneticWithCoil.get_magnetic());
            // painter.export_svg();
        }
    }

    TEST(Test_CoilAdviser_Random_2) {
        srand (time(NULL));

        std::vector<double> turnsRatios;
        int64_t numberStacks = 1;

        std::vector<int64_t> numberTurns = {24, 78, 76};
        std::vector<int64_t> numberParallels = {1, 1, 1};

        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = 507026;
        double magnetizingInductance = 10e-6;
        double temperature = 25;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double peakToPeak = 10;
        double dutyCycle = 0.5;
        double dcCurrent = 0;
        auto coreShapeNames = OpenMagnetics::get_shape_names(false);
        std::string coreShapeName;
        OpenMagnetics::MagneticWrapper magnetic;

        auto gapping = OpenMagneticsTesting::get_residual_gap();
        magnetic = OpenMagneticsTesting::get_quick_magnetic("ETD 19", gapping, numberTurns, numberStacks, "3C91");


        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        OpenMagnetics::MasWrapper masMagnetic;
        inputs.process_waveforms();
        masMagnetic.set_inputs(inputs);
        masMagnetic.set_magnetic(magnetic);

        OpenMagnetics::CoilAdviser coilAdviser;
        auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

        if (masMagneticsWithCoil.size() > 0) {
            auto masMagneticWithCoil = masMagneticsWithCoil[0].first;
            OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
            std::string filePath = __FILE__;
            auto outputFilePath = filePath.substr(0, filePath.rfind("/")).append("/../output/");
            auto outFile = outputFilePath;
            std::string filename = "Test_CoilAdviser" + std::to_string(std::rand()) + ".svg";
            outFile.append(filename);
            // OpenMagnetics::Painter painter(outFile);

            // painter.paint_core(masMagneticWithCoil.get_magnetic());
            // painter.paint_bobbin(masMagneticWithCoil.get_magnetic());
            // painter.paint_coil_turns(masMagneticWithCoil.get_magnetic());
            // painter.export_svg();
        }

    }

    TEST(Test_CoilAdviser_Random_3) {
        srand (time(NULL));

        std::vector<double> turnsRatios;
        int64_t numberStacks = 1;

        std::vector<int64_t> numberTurns = {92, 70, 47};
        std::vector<int64_t> numberParallels = {1, 1, 1};

        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = 313655;
        double magnetizingInductance = 0.002571;
        double temperature = 25;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double peakToPeak = 18;
        double dutyCycle = 0.88;
        double dcCurrent = 6;
        auto coreShapeNames = OpenMagnetics::get_shape_names(false);
        std::string coreShapeName;
        OpenMagnetics::MagneticWrapper magnetic;

        OpenMagnetics::CoilAlignment turnsAlignment = magic_enum::enum_cast<OpenMagnetics::CoilAlignment>(3).value();
        OpenMagnetics::WindingOrientation windingOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(0).value();
        OpenMagnetics::WindingOrientation layersOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(2).value();
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.00756);
        magnetic = OpenMagneticsTesting::get_quick_magnetic("U 81/39/20", gapping, numberTurns, numberStacks, "3C91");

        magnetic.get_mutable_coil().set_turns_alignment(turnsAlignment);
        magnetic.get_mutable_coil().set_winding_orientation(windingOrientation);
        magnetic.get_mutable_coil().set_layers_orientation(layersOrientation);

        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        OpenMagnetics::MasWrapper masMagnetic;
        inputs.process_waveforms();
        masMagnetic.set_inputs(inputs);
        masMagnetic.set_magnetic(magnetic);

        OpenMagnetics::CoilAdviser coilAdviser;
        auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

        if (masMagneticsWithCoil.size() > 0) {
            auto masMagneticWithCoil = masMagneticsWithCoil[0].first;
            OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
            std::string filePath = __FILE__;
            auto outputFilePath = filePath.substr(0, filePath.rfind("/")).append("/../output/");
            auto outFile = outputFilePath;
            std::string filename = "Test_CoilAdviser" + std::to_string(std::rand()) + ".svg";
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(masMagneticWithCoil.get_magnetic());
            painter.paint_bobbin(masMagneticWithCoil.get_magnetic());
            painter.paint_coil_turns(masMagneticWithCoil.get_magnetic());
            painter.export_svg();
        }
    }

    TEST(Test_CoilAdviser_Random_4) {
        srand (time(NULL));

        std::vector<double> turnsRatios;
        int64_t numberStacks = 1;

        std::vector<int64_t> numberTurns = {28};
        std::vector<int64_t> numberParallels = {76};
        std::vector<OpenMagnetics::IsolationSide> isolationSides = {OpenMagnetics::IsolationSide::OCTONARY, OpenMagnetics::IsolationSide::QUINARY};

        double frequency = 837961;
        double magnetizingInductance = 0.007191;
        double temperature = 25;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double peakToPeak = 3;
        double dutyCycle = 0.18;
        double dcCurrent = 28;
        auto coreShapeNames = OpenMagnetics::get_shape_names(false);
        std::string coreShapeName;
        OpenMagnetics::MagneticWrapper magnetic;

        OpenMagnetics::CoilAlignment turnsAlignment = magic_enum::enum_cast<OpenMagnetics::CoilAlignment>(2).value();
        OpenMagnetics::WindingOrientation windingOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(0).value();
        OpenMagnetics::WindingOrientation layersOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(0).value();
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.004025, 3);
        magnetic = OpenMagneticsTesting::get_quick_magnetic("UR 39/35/15", gapping, numberTurns, numberStacks, "3C91");

        magnetic.get_mutable_coil().set_turns_alignment(turnsAlignment);
        magnetic.get_mutable_coil().set_winding_orientation(windingOrientation);
        magnetic.get_mutable_coil().set_layers_orientation(layersOrientation);

        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        OpenMagnetics::MasWrapper masMagnetic;
        inputs.process_waveforms();
        masMagnetic.set_inputs(inputs);
        masMagnetic.set_magnetic(magnetic);

        OpenMagnetics::CoilAdviser coilAdviser;
        auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

        if (masMagneticsWithCoil.size() > 0) {
            auto masMagneticWithCoil = masMagneticsWithCoil[0].first;
            OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
            std::string filePath = __FILE__;
            auto outputFilePath = filePath.substr(0, filePath.rfind("/")).append("/../output/");
            auto outFile = outputFilePath;
            std::string filename = "Test_CoilAdviser" + std::to_string(std::rand()) + ".svg";
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(masMagneticWithCoil.get_magnetic());
            painter.paint_bobbin(masMagneticWithCoil.get_magnetic());
            painter.paint_coil_turns(masMagneticWithCoil.get_magnetic());
            painter.export_svg();
        }
    }
}