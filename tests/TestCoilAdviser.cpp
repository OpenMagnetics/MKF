#include "Insulation.h"
#include "Painter.h"
#include "CoreAdviser.h"
#include "CoilAdviser.h"
#include "InputsWrapper.h"
#include "TestingUtils.h"
#include "Settings.h"

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
        auto withstandVoltageForWires = coilAdviser.get_solid_insulation_requirements_for_wires(inputs, {0, 1}, 1);
        CHECK_EQUAL(2UL, withstandVoltageForWires.size());
        CHECK_EQUAL(2UL, withstandVoltageForWires[0].size());
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
        auto withstandVoltageForWires = coilAdviser.get_solid_insulation_requirements_for_wires(inputs, {0, 1}, 1);
        CHECK_EQUAL(1UL, withstandVoltageForWires.size());
        CHECK_EQUAL(2UL, withstandVoltageForWires[0].size());

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
        auto withstandVoltageForWires = coilAdviser.get_solid_insulation_requirements_for_wires(inputs, {0, 1}, 1);
        CHECK_EQUAL(3UL, withstandVoltageForWires.size());
        CHECK_EQUAL(2UL, withstandVoltageForWires[0].size());

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
        auto withstandVoltageForWires = coilAdviser.get_solid_insulation_requirements_for_wires(inputs, {0, 1}, 1);
        CHECK_EQUAL(3UL, withstandVoltageForWires.size());
        CHECK_EQUAL(2UL, withstandVoltageForWires[0].size());

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
        auto withstandVoltageForWires = coilAdviser.get_solid_insulation_requirements_for_wires(inputs, {0, 1}, 1);
        CHECK_EQUAL(3UL, withstandVoltageForWires.size());
        CHECK_EQUAL(2UL, withstandVoltageForWires[0].size());

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
        auto withstandVoltageForWires = coilAdviser.get_solid_insulation_requirements_for_wires(inputs, {0, 1}, 1);
        CHECK_EQUAL(2UL, withstandVoltageForWires.size());
        CHECK_EQUAL(2UL, withstandVoltageForWires[0].size());
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
        auto withstandVoltageForWires = coilAdviser.get_solid_insulation_requirements_for_wires(inputs, {0, 1}, 1);
        CHECK_EQUAL(1UL, withstandVoltageForWires.size());
        CHECK_EQUAL(2UL, withstandVoltageForWires[0].size());

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
        auto withstandVoltageForWires = coilAdviser.get_solid_insulation_requirements_for_wires(inputs, {0, 1}, 1);
        CHECK_EQUAL(3UL, withstandVoltageForWires.size());
        CHECK_EQUAL(2UL, withstandVoltageForWires[0].size());

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
        auto withstandVoltageForWires = coilAdviser.get_solid_insulation_requirements_for_wires(inputs, {0, 1}, 1);
        CHECK_EQUAL(3UL, withstandVoltageForWires.size());
        CHECK_EQUAL(2UL, withstandVoltageForWires[0].size());

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
        auto withstandVoltageForWires = coilAdviser.get_solid_insulation_requirements_for_wires(inputs, {0, 1}, 1);
        CHECK_EQUAL(3UL, withstandVoltageForWires.size());
        CHECK_EQUAL(2UL, withstandVoltageForWires[0].size());

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
        auto withstandVoltageForWires = coilAdviser.get_solid_insulation_requirements_for_wires(inputs, {0, 1}, 1);
        CHECK_EQUAL(6UL, withstandVoltageForWires.size());
        CHECK_EQUAL(2UL, withstandVoltageForWires[0].size());

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
}

SUITE(CoilAdviser) {
    auto settings = OpenMagnetics::Settings::GetInstance();

    TEST(Test_CoilAdviser_Insulation_No_Margin) {
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.003);
        std::vector<double> turnsRatios;
        int64_t numberStacks = 1;


        double magnetizingInductance = 10e-6;
        double temperature = 25;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double dutyCycle = 0.5;
        double dcCurrent = 0;

        std::vector<int64_t> numberTurns = {82, 5};
        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }
        double frequency = 175590;
        double peakToPeak = 20;

        auto magnetic = OpenMagneticsTesting::get_quick_magnetic("ETD 59", gapping, numberTurns, numberStacks, "3C91");
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        OpenMagnetics::DimensionWithTolerance altitude;
        OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641, OpenMagnetics::InsulationStandards::IEC_623681};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        auto insulationRequirements = OpenMagneticsTesting::get_quick_insulation_requirements(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards);
        inputs.get_mutable_design_requirements().set_insulation(insulationRequirements);

        OpenMagnetics::MasWrapper masMagnetic;
        inputs.process_waveforms();
        masMagnetic.set_inputs(inputs);
        masMagnetic.set_magnetic(magnetic);

        auto settings = OpenMagnetics::Settings::GetInstance();
        settings->set_coil_allow_margin_tape(false);
        settings->set_coil_allow_insulated_wire(true);
        settings->set_coil_try_rewind(false);
        settings->set_coil_adviser_maximum_number_wires(1000);

        OpenMagnetics::CoilAdviser coilAdviser;
        auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 1);


        CHECK(masMagneticsWithCoil.size() > 0);

        int currentIndex = 0;
        for (auto& masMagneticWithCoil : masMagneticsWithCoil) {
            OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
            auto bobbin = masMagneticWithCoil.get_mutable_magnetic().get_bobbin();
            auto aux = bobbin.get_winding_window_dimensions();
            auto bobbinWindingWindowHeight = aux[1];

            auto standardCoordinator = OpenMagnetics::InsulationCoordinator();

            auto creepageDistance = standardCoordinator.calculate_creepage_distance(inputs);
            auto clearance = standardCoordinator.calculate_clearance(inputs);

            auto primarySectionHeight = masMagneticWithCoil.get_mutable_magnetic().get_mutable_coil().get_sections_description().value()[0].get_dimensions()[1];

            int64_t numberMaximumLayers = 0;
            int64_t maximumGrade = 0;
            for (auto wire : masMagneticWithCoil.get_mutable_magnetic().get_wires()) {

                OpenMagnetics::InsulationWireCoating wireCoating;
                if (wire.get_type() == OpenMagnetics::WireType::LITZ) {
                    auto strand = wire.resolve_strand();
                    wireCoating = OpenMagnetics::WireWrapper::resolve_coating(strand).value();
                }
                else {
                    wireCoating = wire.resolve_coating().value();
                }
                if (wireCoating.get_number_layers()) {
                    numberMaximumLayers = std::max(numberMaximumLayers, wireCoating.get_number_layers().value());
                }
                else if (wireCoating.get_grade()) {
                    maximumGrade = std::max(maximumGrade, wireCoating.get_grade().value());
                }
                else {
                    throw std::invalid_argument("Wire must have either layers or grade: " + wire.get_name().value());
                }
            }

            CHECK(numberMaximumLayers >= 3 || maximumGrade >= 3);
            CHECK((bobbinWindingWindowHeight - primarySectionHeight) < std::max(creepageDistance, clearance));

            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_CoilAdviser_No_Margin_" + std::to_string(currentIndex) + ".svg";
            currentIndex++;
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
            painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
            painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
            // painter.paint_coil_sections(masMagneticWithCoil.get_mutable_magnetic());
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_CoilAdviser_Insulation_Margin) {
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.003);
        std::vector<double> turnsRatios;
        int64_t numberStacks = 1;


        double magnetizingInductance = 10e-6;
        double temperature = 25;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double dutyCycle = 0.5;
        double dcCurrent = 0;

        std::vector<int64_t> numberTurns = {82, 55};
        std::vector<int64_t> numberParallels = {1, 1};
        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }
        double frequency = 75590;
        double peakToPeak = 13;

        auto magnetic = OpenMagneticsTesting::get_quick_magnetic("ETD 54", gapping, numberTurns, numberStacks, "3C91");
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        OpenMagnetics::DimensionWithTolerance altitude;
        OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641, OpenMagnetics::InsulationStandards::IEC_623681};
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_IV;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        auto insulationRequirements = OpenMagneticsTesting::get_quick_insulation_requirements(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards);
        inputs.get_mutable_design_requirements().set_insulation(insulationRequirements);

        OpenMagnetics::MasWrapper masMagnetic;
        inputs.process_waveforms();
        masMagnetic.set_inputs(inputs);
        masMagnetic.set_magnetic(magnetic);

        auto settings = OpenMagnetics::Settings::GetInstance();
        settings->set_coil_allow_margin_tape(true);
        settings->set_coil_allow_insulated_wire(false);
        settings->set_coil_try_rewind(false);

        OpenMagnetics::CoilAdviser coilAdviser;
        auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 1);


        CHECK(masMagneticsWithCoil.size() > 0);

        int currentIndex = 0;
        for (auto& masMagneticWithCoil : masMagneticsWithCoil) {
            OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
            auto bobbin = masMagneticWithCoil.get_mutable_magnetic().get_bobbin();
            auto aux = bobbin.get_winding_window_dimensions();
            auto bobbinWindingWindowHeight = aux[1];

            auto standardCoordinator = OpenMagnetics::InsulationCoordinator();

            auto creepageDistance = standardCoordinator.calculate_creepage_distance(inputs);
            auto clearance = standardCoordinator.calculate_clearance(inputs);

            auto primarySectionHeight = masMagneticWithCoil.get_mutable_magnetic().get_mutable_coil().get_sections_description().value()[0].get_dimensions()[1];

            int64_t numberMaximumLayers = 0;
            int64_t maximumGrade = 0;
            for (auto wire : masMagneticWithCoil.get_mutable_magnetic().get_wires()) {

                auto wireCoating = wire.resolve_coating().value();
                if (wireCoating.get_number_layers()) {
                    numberMaximumLayers = std::max(numberMaximumLayers, wireCoating.get_number_layers().value());
                }
                else if (wireCoating.get_grade()) {
                    maximumGrade = std::max(maximumGrade, wireCoating.get_grade().value());
                }
            }

            CHECK(numberMaximumLayers < 3);
            CHECK((bobbinWindingWindowHeight - primarySectionHeight) >= std::max(creepageDistance, clearance));

            CHECK_CLOSE(masMagneticWithCoil.get_mutable_magnetic().get_mutable_coil().get_sections_description().value()[0].get_margin().value()[0], std::max(creepageDistance, clearance) / 2, 0.00001);
            CHECK_CLOSE(masMagneticWithCoil.get_mutable_magnetic().get_mutable_coil().get_sections_description().value()[0].get_margin().value()[1], std::max(creepageDistance, clearance) / 2, 0.00001);

            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_CoilAdviser_Margin_" + std::to_string(currentIndex) + ".svg";
            currentIndex++;
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
            painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
            painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
            // painter.paint_coil_sections(masMagneticWithCoil.get_mutable_magnetic());
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_CoilAdviser_Random) {
        auto settings = OpenMagnetics::Settings::GetInstance();
        settings->reset();
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
                int64_t numberPhysicalTurnsThisWinding = std::rand() % 300 + 1UL;
                int64_t numberTurnsThisWinding = std::rand() % 100 + 1L;
                numberTurns.push_back(numberTurnsThisWinding);
                isolationSides.push_back(OpenMagnetics::get_isolation_side_from_index(static_cast<size_t>(std::rand() % 10 + 1L)));
                int64_t numberParallelsThisWinding = std::max(1.0, std::ceil(double(numberPhysicalTurnsThisWinding) / numberTurnsThisWinding));
                numberParallels.push_back(numberParallelsThisWinding);
                numberPhysicalTurns = std::min(numberPhysicalTurns, numberTurnsThisWinding * numberParallelsThisWinding);
            }
            for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
                turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
            }

            double frequency = std::rand() % 1000000 + 10000L;
            double magnetizingInductance = double(std::rand() % 10000) *  1e-6;
            double temperature = 25;
            OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::SINUSOIDAL;
            double peakToPeak = std::rand() % 30 + 1L;
            double dutyCycle = double(std::rand() % 99 + 1L) / 100;
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
            int windingOrientationIndex = std::rand() % 2;
            int layersOrientationIndex = std::rand() % 2;
            OpenMagnetics::CoilAlignment turnsAlignment = magic_enum::enum_cast<OpenMagnetics::CoilAlignment>(turnsAlignmentIndex).value();
            OpenMagnetics::WindingOrientation windingOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(windingOrientationIndex).value();
            OpenMagnetics::WindingOrientation layersOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(layersOrientationIndex).value();


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
            try {
                auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

                if (masMagneticsWithCoil.size() > 0) {
                    auto masMagneticWithCoil = masMagneticsWithCoil[0];
                    if (!masMagneticWithCoil.get_magnetic().get_coil().get_turns_description()) {
                        continue;
                    }
                    bool result = OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
                    if (!result) {
                        throw std::invalid_argument("Just to trigger prints and return");
                    }
                    count--;
                }
            }
            catch(...) {
                for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
                    std::cout << "numberTurns: " << numberTurns[windingIndex] << std::endl;
                }
                for (size_t windingIndex = 0; windingIndex < numberParallels.size(); ++windingIndex) {
                    std::cout << "numberParallels: " << numberParallels[windingIndex] << std::endl;
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
                CHECK(false);
                return;
            }
        }
        settings->reset();
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
            auto masMagneticWithCoil = masMagneticsWithCoil[0];
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
            auto masMagneticWithCoil = masMagneticsWithCoil[0];
            OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_CoilAdviser" + std::to_string(std::rand()) + ".svg";
            outFile.append(filename);
            // OpenMagnetics::Painter painter(outFile);

            // painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
            // painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
            // painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
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
            auto masMagneticWithCoil = masMagneticsWithCoil[0];
            OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_CoilAdviser" + std::to_string(std::rand()) + ".svg";
            outFile.append(filename);
            // OpenMagnetics::Painter painter(outFile);

            // painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
            // painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
            // painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
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
        OpenMagnetics::WindingOrientation layersOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(1).value();
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
            auto masMagneticWithCoil = masMagneticsWithCoil[0];
            OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_CoilAdviser" + std::to_string(std::rand()) + ".svg";
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
            painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
            painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
            painter.export_svg();
        }
    }

    TEST(Test_CoilAdviser_Random_4) {
        srand (time(NULL));
        auto settings = OpenMagnetics::Settings::GetInstance();
        settings->set_coil_wind_even_if_not_fit(false);

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
            auto masMagneticWithCoil = masMagneticsWithCoil[0];
            OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_CoilAdviser" + std::to_string(std::rand()) + ".svg";
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
            painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
            painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
            painter.export_svg();
        }
        settings->reset();
    }

    TEST(Test_CoilAdviser_Random_5) {
        srand (time(NULL));

        int64_t numberStacks = 1;

        std::vector<int64_t> numberTurns = {15, 36, 87, 60};
        std::vector<int64_t> numberParallels = {1, 1, 1, 1};
        std::vector<double> turnsRatios = {15. / 36, 15. / 87, 15. / 60};
        std::vector<OpenMagnetics::IsolationSide> isolationSides = {OpenMagnetics::IsolationSide::DENARY,
                                                                    OpenMagnetics::IsolationSide::NONARY,
                                                                    OpenMagnetics::IsolationSide::QUATERNARY,
                                                                    OpenMagnetics::IsolationSide::OCTONARY};

        double frequency = 592535;
        double magnetizingInductance = 0.002575;
        double temperature = 23;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double peakToPeak = 3;
        double dutyCycle = 0.61;
        double dcCurrent = 5;
        auto coreShapeNames = OpenMagnetics::get_shape_names(false);
        std::string coreShapeName;
        OpenMagnetics::MagneticWrapper magnetic;

        OpenMagnetics::CoilAlignment turnsAlignment = magic_enum::enum_cast<OpenMagnetics::CoilAlignment>(3).value();
        OpenMagnetics::WindingOrientation windingOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(0).value();
        OpenMagnetics::WindingOrientation layersOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(1).value();
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.003175, 1);
        magnetic = OpenMagneticsTesting::get_quick_magnetic("E 114/46/26", gapping, numberTurns, numberStacks, "3C91");

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
            auto masMagneticWithCoil = masMagneticsWithCoil[0];
            OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_CoilAdviser" + std::to_string(std::rand()) + ".svg";
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
            painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
            painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
            painter.export_svg();
        }
    }

    TEST(Test_CoilAdviser_Random_6) {
        srand (time(NULL));

        int64_t numberStacks = 1;

        std::vector<int64_t> numberTurns = {11, 82};
        std::vector<int64_t> numberParallels = {3, 3};
        std::vector<double> turnsRatios = {11. / 82};
        std::vector<OpenMagnetics::IsolationSide> isolationSides = {OpenMagnetics::IsolationSide::DENARY,
                                                                    OpenMagnetics::IsolationSide::NONARY};

        double frequency = 617645;
        double magnetizingInductance = 0.009088;
        double temperature = 23;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double peakToPeak = 26;
        double dutyCycle = 0.38;
        double dcCurrent = 16;
        auto coreShapeNames = OpenMagnetics::get_shape_names(false);
        std::string coreShapeName;
        OpenMagnetics::MagneticWrapper magnetic;

        OpenMagnetics::CoilAlignment turnsAlignment = magic_enum::enum_cast<OpenMagnetics::CoilAlignment>(1).value();
        OpenMagnetics::WindingOrientation windingOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(0).value();
        OpenMagnetics::WindingOrientation layersOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(1).value();
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.001227);
        magnetic = OpenMagneticsTesting::get_quick_magnetic("E 72/28/19", gapping, numberTurns, numberStacks, "3C91");

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
            auto masMagneticWithCoil = masMagneticsWithCoil[0];
            OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_CoilAdviser" + std::to_string(std::rand()) + ".svg";
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
            painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
            painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
            painter.export_svg();
        }
    }

    TEST(Test_CoilAdviser_Random_7) {
        srand (time(NULL));

        int64_t numberStacks = 1;

        std::vector<int64_t> numberTurns = {60, 18, 10};
        std::vector<int64_t> numberParallels = {1, 3, 3};
        std::vector<double> turnsRatios = {60. / 18, 60. / 10};
        std::vector<OpenMagnetics::IsolationSide> isolationSides = {OpenMagnetics::IsolationSide::SECONDARY,
                                                                    OpenMagnetics::IsolationSide::NONARY,
                                                                    OpenMagnetics::IsolationSide::TERTIARY};

        double frequency = 95989;
        double magnetizingInductance = 0.009266;
        double temperature = 23;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double peakToPeak = 9;
        double dutyCycle = 0.57;
        double dcCurrent = 1;
        auto coreShapeNames = OpenMagnetics::get_shape_names(false);
        std::string coreShapeName;
        OpenMagnetics::MagneticWrapper magnetic;

        OpenMagnetics::CoilAlignment turnsAlignment = magic_enum::enum_cast<OpenMagnetics::CoilAlignment>(1).value();
        OpenMagnetics::WindingOrientation windingOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(1).value();
        OpenMagnetics::WindingOrientation layersOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(0).value();
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.007023);
        magnetic = OpenMagneticsTesting::get_quick_magnetic("E 35", gapping, numberTurns, numberStacks, "3C91");

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
            auto masMagneticWithCoil = masMagneticsWithCoil[0];
            OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_CoilAdviser" + std::to_string(std::rand()) + ".svg";
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
            painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
            painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
            painter.export_svg();
        }
    }

    TEST(Test_CoilAdviser_Random_8) {
        srand (time(NULL));

        int64_t numberStacks = 1;

        std::vector<int64_t> numberTurns = {36, 3, 8};
        std::vector<int64_t> numberParallels = {2, 1, 2};
        std::vector<double> turnsRatios = {36. / 3, 36. / 8};
        std::vector<OpenMagnetics::IsolationSide> isolationSides = {OpenMagnetics::IsolationSide::UNDENARY,
                                                                    OpenMagnetics::IsolationSide::OCTONARY,
                                                                    OpenMagnetics::IsolationSide::DENARY};

        double frequency = 632226;
        double magnetizingInductance = 0.001529;
        double temperature = 23;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double peakToPeak = 3;
        double dutyCycle = 0.68;
        double dcCurrent = 7;
        auto coreShapeNames = OpenMagnetics::get_shape_names(false);
        std::string coreShapeName;
        OpenMagnetics::MagneticWrapper magnetic;

        OpenMagnetics::CoilAlignment turnsAlignment = magic_enum::enum_cast<OpenMagnetics::CoilAlignment>(0).value();
        OpenMagnetics::WindingOrientation windingOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(0).value();
        OpenMagnetics::WindingOrientation layersOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(0).value();
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.009828);
        magnetic = OpenMagneticsTesting::get_quick_magnetic("E 80/38/20", gapping, numberTurns, numberStacks, "3C91");

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
            auto masMagneticWithCoil = masMagneticsWithCoil[0];
            OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_CoilAdviser" + std::to_string(std::rand()) + ".svg";
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
            painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
            painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
            painter.export_svg();
        }
    }

    TEST(Test_CoilAdviser_Random_9) {
        srand (time(NULL));

        int64_t numberStacks = 1;

        std::vector<int64_t> numberTurns = {36, 55, 96};
        std::vector<int64_t> numberParallels = {2, 1, 2};
        std::vector<double> turnsRatios = {36. / 55, 36. / 96};
        std::vector<OpenMagnetics::IsolationSide> isolationSides = {OpenMagnetics::IsolationSide::SENARY,
                                                                    OpenMagnetics::IsolationSide::NONARY,
                                                                    OpenMagnetics::IsolationSide::NONARY};

        double frequency = 632226;
        double magnetizingInductance = 0.001529;
        double temperature = 23;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double peakToPeak = 3;
        double dutyCycle = 0.68;
        double dcCurrent = 7;
        auto coreShapeNames = OpenMagnetics::get_shape_names(false);
        std::string coreShapeName;
        OpenMagnetics::MagneticWrapper magnetic;

        OpenMagnetics::CoilAlignment turnsAlignment = magic_enum::enum_cast<OpenMagnetics::CoilAlignment>(2).value();
        OpenMagnetics::WindingOrientation windingOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(0).value();
        OpenMagnetics::WindingOrientation layersOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(0).value();
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.009828);
        magnetic = OpenMagneticsTesting::get_quick_magnetic("U 79/129/31", gapping, numberTurns, numberStacks, "3C91");

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
            auto masMagneticWithCoil = masMagneticsWithCoil[0];
            OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_CoilAdviser" + std::to_string(std::rand()) + ".svg";
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
            painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
            painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
            painter.export_svg();
        }
    }

    TEST(Test_CoilAdviser_Random_10) {
        srand (time(NULL));

        int64_t numberStacks = 1;

        std::vector<int64_t> numberTurns = {49, 80, 78, 1};
        std::vector<int64_t> numberParallels = {3, 1, 1, 2};
        std::vector<double> turnsRatios = {49. / 80, 49. / 78, 49};
        std::vector<OpenMagnetics::IsolationSide> isolationSides = {OpenMagnetics::IsolationSide::TERTIARY,
                                                                    OpenMagnetics::IsolationSide::SENARY,
                                                                    OpenMagnetics::IsolationSide::DENARY,
                                                                    OpenMagnetics::IsolationSide::NONARY};

        double frequency = 660462;
        double magnetizingInductance = 0.006606;
        double temperature = 23;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double peakToPeak = 24;
        double dutyCycle = 0.28;
        double dcCurrent = 2;
        auto coreShapeNames = OpenMagnetics::get_shape_names(false);
        std::string coreShapeName;
        OpenMagnetics::MagneticWrapper magnetic;

        OpenMagnetics::CoilAlignment turnsAlignment = magic_enum::enum_cast<OpenMagnetics::CoilAlignment>(2).value();
        OpenMagnetics::WindingOrientation windingOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(1).value();
        OpenMagnetics::WindingOrientation layersOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(1).value();
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.002048, 3);
        magnetic = OpenMagneticsTesting::get_quick_magnetic("E 50/15", gapping, numberTurns, numberStacks, "3C91");

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
            auto masMagneticWithCoil = masMagneticsWithCoil[0];
            // OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_CoilAdviser" + std::to_string(std::rand()) + ".svg";
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
            painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
            painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
            painter.export_svg();
        }
    }

    TEST(Test_CoilAdviser_Random_11) {
        srand (time(NULL));

        int64_t numberStacks = 1;

        std::vector<int64_t> numberTurns = {72};
        std::vector<int64_t> numberParallels = {4};
        std::vector<double> turnsRatios = {};
        std::vector<OpenMagnetics::IsolationSide> isolationSides = {OpenMagnetics::IsolationSide::NONARY};

        double frequency = 821021;
        double magnetizingInductance = 8.6e-05;
        double temperature = 23;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double peakToPeak = 19;
        double dutyCycle = 0.14;
        double dcCurrent = 0;
        auto coreShapeNames = OpenMagnetics::get_shape_names(false);
        std::string coreShapeName;
        OpenMagnetics::MagneticWrapper magnetic;

        OpenMagnetics::CoilAlignment turnsAlignment = magic_enum::enum_cast<OpenMagnetics::CoilAlignment>(3).value();
        OpenMagnetics::WindingOrientation windingOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(0).value();
        OpenMagnetics::WindingOrientation layersOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(1).value();
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.001424, 1);
        magnetic = OpenMagneticsTesting::get_quick_magnetic("U 15/11/6", gapping, numberTurns, numberStacks, "3C91");

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
            auto masMagneticWithCoil = masMagneticsWithCoil[0];
            // OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_CoilAdviser" + std::to_string(std::rand()) + ".svg";
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
            painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
            painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
            painter.export_svg();
        }
    }

    TEST(Test_CoilAdviser_Random_12) {
        srand (time(NULL));

        int64_t numberStacks = 1;

        std::vector<int64_t> numberTurns = {53, 100, 80, 98};
        std::vector<int64_t> numberParallels = {5, 2, 3, 2};
        std::vector<double> turnsRatios = {53.0 / 100, 53.0 / 80, 53.0 / 98};
        std::vector<OpenMagnetics::IsolationSide> isolationSides = {OpenMagnetics::IsolationSide::OCTONARY,
                                                                    OpenMagnetics::IsolationSide::SENARY,
                                                                    OpenMagnetics::IsolationSide::QUINARY,
                                                                    OpenMagnetics::IsolationSide::SENARY};


        double frequency = 460425;
        double magnetizingInductance = 0.005275;
        double temperature = 23;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double peakToPeak = 28;
        double dutyCycle = 0.73;
        double dcCurrent = 5;
        auto coreShapeNames = OpenMagnetics::get_shape_names(false);
        std::string coreShapeName;
        OpenMagnetics::MagneticWrapper magnetic;

        OpenMagnetics::CoilAlignment turnsAlignment = magic_enum::enum_cast<OpenMagnetics::CoilAlignment>(0).value();
        OpenMagnetics::WindingOrientation windingOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(1).value();
        OpenMagnetics::WindingOrientation layersOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(0).value();

        auto gapping = OpenMagneticsTesting::get_spacer_gap(0.006456);
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

        inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        OpenMagnetics::MasWrapper masMagnetic;
        inputs.process_waveforms();
        masMagnetic.set_inputs(inputs);
        masMagnetic.set_magnetic(magnetic);

        OpenMagnetics::CoilAdviser coilAdviser;
        auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

        if (masMagneticsWithCoil.size() > 0) {
            auto masMagneticWithCoil = masMagneticsWithCoil[0];
            OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            std::string filename = "Test_CoilAdviser_Random_12.svg";
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
            painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
            painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
            painter.export_svg();
        }
    }
}