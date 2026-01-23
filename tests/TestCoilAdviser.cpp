#include "RandomUtils.h"
#include <source_location>
#include "constructive_models/Insulation.h"
#include "constructive_models/CorePiece.h"
#include "support/Painter.h"
#include "advisers/CoreAdviser.h"
#include "advisers/CoilAdviser.h"
#include "processors/Inputs.h"
#include <magic_enum.hpp>
#include "TestingUtils.h"
#include "support/Settings.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;

namespace {

auto standardCoordinator = InsulationCoordinator();
auto overvoltageCategory = OvervoltageCategory::OVC_II;
auto cti = Cti::GROUP_I;
double maximumVoltageRms = 666;
double maximumVoltagePeak = 800;
double frequency = 30000;
DimensionWithTolerance altitude;
DimensionWithTolerance mainSupplyVoltage;
auto pollutionDegree = PollutionDegree::P1;


TEST_CASE("Test_CoilAdviser_Wires_Withstand_Voltage_Same_Isolation_Sides_Basic_No_FIW", "[adviser][coil-adviser][insulation][smoke-test]") {

    auto standards = std::vector<InsulationStandards>{InsulationStandards::IEC_606641};
    altitude.set_maximum(2000);
    mainSupplyVoltage.set_nominal(400);
    auto insulationType = InsulationType::BASIC;
    OpenMagnetics::Inputs inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, WiringTechnology::WOUND);
    DimensionWithTolerance dimensionWithTolerance;
    dimensionWithTolerance.set_nominal(1);
    inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
    inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<IsolationSide>{IsolationSide::PRIMARY, IsolationSide::PRIMARY});

    CoilAdviser coilAdviser;
    auto withstandVoltageForWires = InsulationCoordinator::get_solid_insulation_requirements_for_wires(inputs, {0, 1}, 1);
    REQUIRE(2UL == withstandVoltageForWires.size());
    REQUIRE(2UL == withstandVoltageForWires[0].size());
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_grade());

    REQUIRE(withstandVoltageForWires[0][0].get_minimum_grade().value() == 1);
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_breakdown_voltage() > 0);
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_breakdown_voltage() > 0);
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_number_layers().value() == 1);
    REQUIRE(!withstandVoltageForWires[1][0].get_minimum_grade());
}

TEST_CASE("Test_CoilAdviser_Wires_Withstand_Voltage_Different_Isolation_Sides_Functional_No_FIW", "[adviser][coil-adviser][insulation][smoke-test]") {
    auto standards = std::vector<InsulationStandards>{InsulationStandards::IEC_606641};
    altitude.set_maximum(2000);
    mainSupplyVoltage.set_nominal(400);
    auto insulationType = InsulationType::FUNCTIONAL;
    OpenMagnetics::Inputs inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, WiringTechnology::WOUND);
    DimensionWithTolerance dimensionWithTolerance;
    dimensionWithTolerance.set_nominal(1);
    inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
    inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<IsolationSide>{IsolationSide::PRIMARY, IsolationSide::SECONDARY});

    CoilAdviser coilAdviser;
    auto withstandVoltageForWires = InsulationCoordinator::get_solid_insulation_requirements_for_wires(inputs, {0, 1}, 1);
    REQUIRE(1UL == withstandVoltageForWires.size());
    REQUIRE(2UL == withstandVoltageForWires[0].size());

    REQUIRE(withstandVoltageForWires[0][0].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_grade());
}

TEST_CASE("Test_CoilAdviser_Wires_Withstand_Voltage_Different_Isolation_Sides_Basic_No_FIW", "[adviser][coil-adviser][insulation][smoke-test]") {
    auto standards = std::vector<InsulationStandards>{InsulationStandards::IEC_606641};
    altitude.set_maximum(2000);
    mainSupplyVoltage.set_nominal(400);
    auto insulationType = InsulationType::BASIC;
    OpenMagnetics::Inputs inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, WiringTechnology::WOUND);
    DimensionWithTolerance dimensionWithTolerance;
    dimensionWithTolerance.set_nominal(1);
    inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
    inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<IsolationSide>{IsolationSide::PRIMARY, IsolationSide::SECONDARY});

    CoilAdviser coilAdviser;
    auto withstandVoltageForWires = InsulationCoordinator::get_solid_insulation_requirements_for_wires(inputs, {0, 1}, 1);
    REQUIRE(3UL == withstandVoltageForWires.size());
    REQUIRE(2UL == withstandVoltageForWires[0].size());

    REQUIRE(withstandVoltageForWires[0][0].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_grade());

    REQUIRE(withstandVoltageForWires[1][0].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_breakdown_voltage() == 6000);
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_grade());
    REQUIRE(!withstandVoltageForWires[1][1].get_minimum_grade());

    REQUIRE(withstandVoltageForWires[2][0].get_minimum_breakdown_voltage() == 6000);
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[2][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[2][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_number_layers().value() == 1);
    REQUIRE(!withstandVoltageForWires[2][0].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_grade());
}

TEST_CASE("Test_CoilAdviser_Wires_Withstand_Voltage_Different_Isolation_Sides_Supplementary_No_FIW", "[adviser][coil-adviser][insulation][smoke-test]") {
    auto standards = std::vector<InsulationStandards>{InsulationStandards::IEC_606641};
    altitude.set_maximum(2000);
    mainSupplyVoltage.set_nominal(400);
    auto insulationType = InsulationType::SUPPLEMENTARY;
    OpenMagnetics::Inputs inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, WiringTechnology::WOUND);
    DimensionWithTolerance dimensionWithTolerance;
    dimensionWithTolerance.set_nominal(1);
    inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
    inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<IsolationSide>{IsolationSide::PRIMARY, IsolationSide::SECONDARY});

    CoilAdviser coilAdviser;
    auto withstandVoltageForWires = InsulationCoordinator::get_solid_insulation_requirements_for_wires(inputs, {0, 1}, 1);
    REQUIRE(3UL == withstandVoltageForWires.size());
    REQUIRE(2UL == withstandVoltageForWires[0].size());

    REQUIRE(withstandVoltageForWires[0][0].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_grade());

    REQUIRE(withstandVoltageForWires[1][0].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_breakdown_voltage() == 6000);
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_grade());
    REQUIRE(!withstandVoltageForWires[1][1].get_minimum_grade());

    REQUIRE(withstandVoltageForWires[2][0].get_minimum_breakdown_voltage() == 6000);
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[2][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[2][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_number_layers().value() == 1);
    REQUIRE(!withstandVoltageForWires[2][0].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_grade());
}

TEST_CASE("Test_CoilAdviser_Wires_Withstand_Voltage_Different_Isolation_Sides_Reinforced_No_FIW", "[adviser][coil-adviser][insulation][smoke-test]") {
    auto standards = std::vector<InsulationStandards>{InsulationStandards::IEC_606641};
    altitude.set_maximum(2000);
    mainSupplyVoltage.set_nominal(400);
    auto insulationType = InsulationType::REINFORCED;
    OpenMagnetics::Inputs inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, WiringTechnology::WOUND);
    DimensionWithTolerance dimensionWithTolerance;
    dimensionWithTolerance.set_nominal(1);
    inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
    inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<IsolationSide>{IsolationSide::PRIMARY, IsolationSide::SECONDARY});

    CoilAdviser coilAdviser;
    auto withstandVoltageForWires = InsulationCoordinator::get_solid_insulation_requirements_for_wires(inputs, {0, 1}, 1);
    REQUIRE(3UL == withstandVoltageForWires.size());
    REQUIRE(2UL == withstandVoltageForWires[0].size());

    REQUIRE(withstandVoltageForWires[0][0].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_grade());

    REQUIRE(withstandVoltageForWires[1][0].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_breakdown_voltage() == 8000);
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_number_layers().value() == 3);
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_grade());
    REQUIRE(!withstandVoltageForWires[1][1].get_minimum_grade());

    REQUIRE(withstandVoltageForWires[2][0].get_minimum_breakdown_voltage() == 8000);
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[2][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[2][0].get_minimum_number_layers().value() == 3);
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_number_layers().value() == 1);
    REQUIRE(!withstandVoltageForWires[2][0].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_grade());
}

TEST_CASE("Test_CoilAdviser_Wires_Withstand_Voltage_Same_Isolation_Sides_Basic_FIW", "[adviser][coil-adviser][insulation][smoke-test]") {

    auto standards = std::vector<InsulationStandards>{InsulationStandards::IEC_623681};
    altitude.set_maximum(2000);
    mainSupplyVoltage.set_nominal(400);
    auto insulationType = InsulationType::BASIC;
    OpenMagnetics::Inputs inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, WiringTechnology::WOUND);
    DimensionWithTolerance dimensionWithTolerance;
    dimensionWithTolerance.set_nominal(1);
    inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
    inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<IsolationSide>{IsolationSide::PRIMARY, IsolationSide::PRIMARY});

    CoilAdviser coilAdviser;
    auto withstandVoltageForWires = InsulationCoordinator::get_solid_insulation_requirements_for_wires(inputs, {0, 1}, 1);
    REQUIRE(2UL == withstandVoltageForWires.size());
    REQUIRE(2UL == withstandVoltageForWires[0].size());
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_grade());

    REQUIRE(withstandVoltageForWires[0][0].get_minimum_grade().value() == 1);
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_breakdown_voltage() > 0);
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_breakdown_voltage() > 0);
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_grade().value() == 3);
}

TEST_CASE("Test_CoilAdviser_Wires_Withstand_Voltage_Different_Isolation_Sides_Functional_FIW", "[adviser][coil-adviser][insulation][smoke-test]") {
    auto standards = std::vector<InsulationStandards>{InsulationStandards::IEC_623681};
    altitude.set_maximum(2000);
    mainSupplyVoltage.set_nominal(400);
    auto insulationType = InsulationType::FUNCTIONAL;
    OpenMagnetics::Inputs inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, WiringTechnology::WOUND);
    DimensionWithTolerance dimensionWithTolerance;
    dimensionWithTolerance.set_nominal(1);
    inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
    inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<IsolationSide>{IsolationSide::PRIMARY, IsolationSide::SECONDARY});

    CoilAdviser coilAdviser;
    auto withstandVoltageForWires = InsulationCoordinator::get_solid_insulation_requirements_for_wires(inputs, {0, 1}, 1);
    REQUIRE(1UL == withstandVoltageForWires.size());
    REQUIRE(2UL == withstandVoltageForWires[0].size());

    REQUIRE(withstandVoltageForWires[0][0].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_grade());
}

TEST_CASE("Test_CoilAdviser_Wires_Withstand_Voltage_Different_Isolation_Sides_Basic_FIW", "[adviser][coil-adviser][insulation][smoke-test]") {
    auto standards = std::vector<InsulationStandards>{InsulationStandards::IEC_623681};
    altitude.set_maximum(2000);
    mainSupplyVoltage.set_nominal(400);
    auto insulationType = InsulationType::BASIC;
    OpenMagnetics::Inputs inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, WiringTechnology::WOUND);
    DimensionWithTolerance dimensionWithTolerance;
    dimensionWithTolerance.set_nominal(1);
    inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
    inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<IsolationSide>{IsolationSide::PRIMARY, IsolationSide::SECONDARY});

    CoilAdviser coilAdviser;
    auto withstandVoltageForWires = InsulationCoordinator::get_solid_insulation_requirements_for_wires(inputs, {0, 1}, 1);
    REQUIRE(3UL == withstandVoltageForWires.size());
    REQUIRE(2UL == withstandVoltageForWires[0].size());

    REQUIRE(withstandVoltageForWires[0][0].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_grade());

    REQUIRE(withstandVoltageForWires[1][0].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_breakdown_voltage() == 2500);
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_grade().value() == 3);

    REQUIRE(withstandVoltageForWires[2][0].get_minimum_breakdown_voltage() == 2500);
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[2][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[2][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[2][0].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[2][0].get_minimum_grade().value() == 3);
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_grade());
}

TEST_CASE("Test_CoilAdviser_Wires_Withstand_Voltage_Different_Isolation_Sides_Supplementary_FIW", "[adviser][coil-adviser][insulation][smoke-test]") {
    auto standards = std::vector<InsulationStandards>{InsulationStandards::IEC_623681};
    altitude.set_maximum(2000);
    mainSupplyVoltage.set_nominal(400);
    auto insulationType = InsulationType::SUPPLEMENTARY;
    OpenMagnetics::Inputs inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, WiringTechnology::WOUND);
    DimensionWithTolerance dimensionWithTolerance;
    dimensionWithTolerance.set_nominal(1);
    inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
    inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<IsolationSide>{IsolationSide::PRIMARY, IsolationSide::SECONDARY});

    CoilAdviser coilAdviser;
    auto withstandVoltageForWires = InsulationCoordinator::get_solid_insulation_requirements_for_wires(inputs, {0, 1}, 1);
    REQUIRE(3UL == withstandVoltageForWires.size());
    REQUIRE(2UL == withstandVoltageForWires[0].size());

    REQUIRE(withstandVoltageForWires[0][0].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_grade());

    REQUIRE(withstandVoltageForWires[1][0].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_breakdown_voltage() == 2500);
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_grade().value() == 3);

    REQUIRE(withstandVoltageForWires[2][0].get_minimum_breakdown_voltage() == 2500);
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[2][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[2][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[2][0].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[2][0].get_minimum_grade().value() == 3);
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_grade());
}

TEST_CASE("Test_CoilAdviser_Wires_Withstand_Voltage_Different_Isolation_Sides_Reinforced_FIW", "[adviser][coil-adviser][insulation][smoke-test]") {
    auto standards = std::vector<InsulationStandards>{InsulationStandards::IEC_623681};
    altitude.set_maximum(2000);
    mainSupplyVoltage.set_nominal(400);
    auto insulationType = InsulationType::REINFORCED;
    OpenMagnetics::Inputs inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, WiringTechnology::WOUND);
    DimensionWithTolerance dimensionWithTolerance;
    dimensionWithTolerance.set_nominal(1);
    inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
    inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<IsolationSide>{IsolationSide::PRIMARY, IsolationSide::SECONDARY});

    CoilAdviser coilAdviser;
    auto withstandVoltageForWires = InsulationCoordinator::get_solid_insulation_requirements_for_wires(inputs, {0, 1}, 1);
    REQUIRE(3UL == withstandVoltageForWires.size());
    REQUIRE(2UL == withstandVoltageForWires[0].size());

    REQUIRE(withstandVoltageForWires[0][0].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_grade()); 

    REQUIRE(withstandVoltageForWires[1][0].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_breakdown_voltage() == 5000);
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_number_layers().value() == 3);
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_grade().value() == 3);

    REQUIRE(withstandVoltageForWires[2][0].get_minimum_breakdown_voltage() == 5000);
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[2][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[2][0].get_minimum_number_layers().value() == 3);
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[2][0].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[2][0].get_minimum_grade().value() == 3);
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_grade());
}

TEST_CASE("Test_CoilAdviser_Wires_Withstand_Voltage_Different_Isolation_Sides_Double_FIW", "[adviser][coil-adviser][insulation][smoke-test]") {
    auto standards = std::vector<InsulationStandards>{InsulationStandards::IEC_623681};
    altitude.set_maximum(2000);
    mainSupplyVoltage.set_nominal(400);
    auto insulationType = InsulationType::DOUBLE;
    OpenMagnetics::Inputs inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, WiringTechnology::WOUND);
    DimensionWithTolerance dimensionWithTolerance;
    dimensionWithTolerance.set_nominal(1);
    inputs.get_mutable_design_requirements().set_turns_ratios({dimensionWithTolerance});
    inputs.get_mutable_design_requirements().set_isolation_sides(std::vector<IsolationSide>{IsolationSide::PRIMARY, IsolationSide::SECONDARY});

    CoilAdviser coilAdviser;
    auto withstandVoltageForWires = InsulationCoordinator::get_solid_insulation_requirements_for_wires(inputs, {0, 1}, 1);
    REQUIRE(6UL == withstandVoltageForWires.size());
    REQUIRE(2UL == withstandVoltageForWires[0].size());

    REQUIRE(withstandVoltageForWires[0][0].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[0][0].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[0][1].get_minimum_grade());

    REQUIRE(withstandVoltageForWires[1][0].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_breakdown_voltage() == 5000);
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_number_layers().value() == 3);
    REQUIRE(withstandVoltageForWires[1][0].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[1][1].get_minimum_grade().value() == 3);

    REQUIRE(withstandVoltageForWires[2][0].get_minimum_breakdown_voltage() == 5000);
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[2][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[2][0].get_minimum_number_layers().value() == 3);
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[2][0].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[2][0].get_minimum_grade().value() == 3);
    REQUIRE(withstandVoltageForWires[2][1].get_minimum_grade());

    REQUIRE(withstandVoltageForWires[3][0].get_minimum_breakdown_voltage() == 2500);
    REQUIRE(withstandVoltageForWires[3][1].get_minimum_breakdown_voltage() == 2500);
    REQUIRE(withstandVoltageForWires[3][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[3][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[3][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[3][1].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[3][0].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[3][0].get_minimum_grade().value() == 3);
    REQUIRE(withstandVoltageForWires[3][1].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[3][1].get_minimum_grade().value() == 3);

    REQUIRE(withstandVoltageForWires[4][0].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[4][1].get_minimum_breakdown_voltage() == 2500);
    REQUIRE(withstandVoltageForWires[4][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[4][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[4][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[4][1].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[4][0].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[4][1].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[4][1].get_minimum_grade().value() == 3);

    REQUIRE(withstandVoltageForWires[5][0].get_minimum_breakdown_voltage() == 2500);
    REQUIRE(withstandVoltageForWires[5][1].get_minimum_breakdown_voltage() == 0);
    REQUIRE(withstandVoltageForWires[5][0].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[5][0].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[5][1].get_minimum_number_layers());
    REQUIRE(withstandVoltageForWires[5][1].get_minimum_number_layers().value() == 1);
    REQUIRE(withstandVoltageForWires[5][0].get_minimum_grade());
    REQUIRE(withstandVoltageForWires[5][0].get_minimum_grade().value() == 3);
    REQUIRE(withstandVoltageForWires[5][1].get_minimum_grade());
}

TEST_CASE("Test_CoilAdviser_Insulation_No_Margin", "[adviser][coil-adviser][smoke-test]") {
    settings.reset();
    auto gapping = OpenMagneticsTesting::get_ground_gap(0.003);
    std::vector<double> turnsRatios;
    int64_t numberStacks = 1;


    double magnetizingInductance = 10e-6;
    double temperature = 25;
    WaveformLabel waveShape = WaveformLabel::SINUSOIDAL;
    double dutyCycle = 0.5;
    double dcCurrent = 0;

    std::vector<int64_t> numberTurns = {82, 5};
    for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
        turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
    }
    double frequency = 175590;
    double peakToPeak = 20;

    auto magnetic = OpenMagneticsTesting::get_quick_magnetic("ETD 59", gapping, numberTurns, numberStacks, "3C91");
    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                          magnetizingInductance,
                                                                                          temperature,
                                                                                          waveShape,
                                                                                          peakToPeak,
                                                                                          dutyCycle,
                                                                                          dcCurrent,
                                                                                          turnsRatios);

    DimensionWithTolerance altitude;
    DimensionWithTolerance mainSupplyVoltage;
    auto standards = std::vector<InsulationStandards>{InsulationStandards::IEC_606641, InsulationStandards::IEC_623681};
    altitude.set_maximum(2000);
    mainSupplyVoltage.set_nominal(400);
    auto cti = Cti::GROUP_I;
    auto overvoltageCategory = OvervoltageCategory::OVC_IV;
    auto insulationType = InsulationType::BASIC;
    auto pollutionDegree = PollutionDegree::P1;
    auto insulationRequirements = OpenMagneticsTesting::get_quick_insulation_requirements(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards);
    inputs.get_mutable_design_requirements().set_insulation(insulationRequirements);

    OpenMagnetics::Mas masMagnetic;
    inputs.process();
    masMagnetic.set_inputs(inputs);
    masMagnetic.set_magnetic(magnetic);

    settings.set_coil_allow_margin_tape(false);
    settings.set_coil_allow_insulated_wire(true);
    settings.set_coil_try_rewind(false);
    settings.set_coil_adviser_maximum_number_wires(1000);

    CoilAdviser coilAdviser;
    auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 1);


    REQUIRE(masMagneticsWithCoil.size() > 0);

    int currentIndex = 0;
    for (auto& masMagneticWithCoil : masMagneticsWithCoil) {
        OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
        auto bobbin = masMagneticWithCoil.get_mutable_magnetic().get_bobbin();
        auto aux = bobbin.get_winding_window_dimensions();
        auto bobbinWindingWindowHeight = aux[1];

        auto standardCoordinator = InsulationCoordinator();

        auto creepageDistance = standardCoordinator.calculate_creepage_distance(inputs);
        auto clearance = standardCoordinator.calculate_clearance(inputs);

        auto primarySectionHeight = masMagneticWithCoil.get_mutable_magnetic().get_mutable_coil().get_sections_description().value()[0].get_dimensions()[1];

        int64_t numberMaximumLayers = 0;
        int64_t maximumGrade = 0;
        for (auto wire : masMagneticWithCoil.get_mutable_magnetic().get_wires()) {

            InsulationWireCoating wireCoating;
            if (wire.get_type() == WireType::LITZ) {
                auto strand = wire.resolve_strand();
                wireCoating = OpenMagnetics::Wire::resolve_coating(strand).value();
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

        REQUIRE((numberMaximumLayers >= 3 || maximumGrade >= 3));
        REQUIRE((bobbinWindingWindowHeight - primarySectionHeight) < std::max(creepageDistance, clearance));

        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        std::string filename = "Test_CoilAdviser_No_Margin_" + std::to_string(currentIndex) + ".svg";
        currentIndex++;
        outFile.append(filename);
        Painter painter(outFile);

        painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
        // painter.paint_coil_sections(masMagneticWithCoil.get_mutable_magnetic());
        painter.export_svg();
    }
    settings.reset();
}

TEST_CASE("Test_CoilAdviser_Insulation_Margin", "[adviser][coil-adviser][margin][smoke-test]") {
    // OpenMagnetics::set_log_verbosity(2);
    auto gapping = OpenMagneticsTesting::get_ground_gap(0.003);
    std::vector<double> turnsRatios;
    int64_t numberStacks = 1;


    double magnetizingInductance = 10e-6;
    double temperature = 25;
    WaveformLabel waveShape = WaveformLabel::SINUSOIDAL;
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
    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                          magnetizingInductance,
                                                                                          temperature,
                                                                                          waveShape,
                                                                                          peakToPeak,
                                                                                          dutyCycle,
                                                                                          dcCurrent,
                                                                                          turnsRatios);

    DimensionWithTolerance altitude;
    DimensionWithTolerance mainSupplyVoltage;
    auto standards = std::vector<InsulationStandards>{InsulationStandards::IEC_606641, InsulationStandards::IEC_623681};
    altitude.set_maximum(2000);
    mainSupplyVoltage.set_nominal(400);
    auto cti = Cti::GROUP_I;
    auto overvoltageCategory = OvervoltageCategory::OVC_IV;
    auto insulationType = InsulationType::BASIC;
    auto pollutionDegree = PollutionDegree::P1;
    auto insulationRequirements = OpenMagneticsTesting::get_quick_insulation_requirements(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards);
    inputs.get_mutable_design_requirements().set_insulation(insulationRequirements);

    OpenMagnetics::Mas masMagnetic;
    inputs.process();
    masMagnetic.set_inputs(inputs);
    masMagnetic.set_magnetic(magnetic);

    settings.set_coil_allow_margin_tape(true);
    settings.set_coil_allow_insulated_wire(false);
    settings.set_coil_try_rewind(false);

    CoilAdviser coilAdviser;
    auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 1);


    REQUIRE(masMagneticsWithCoil.size() > 0);

    int currentIndex = 0;
    for (auto& masMagneticWithCoil : masMagneticsWithCoil) {
        OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
        auto bobbin = masMagneticWithCoil.get_mutable_magnetic().get_bobbin();
        auto aux = bobbin.get_winding_window_dimensions();
        auto bobbinWindingWindowHeight = aux[1];

        auto standardCoordinator = InsulationCoordinator();

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

        REQUIRE(numberMaximumLayers < 3);
        REQUIRE((bobbinWindingWindowHeight - primarySectionHeight) >= std::max(creepageDistance, clearance));

        REQUIRE_THAT(OpenMagnetics::Coil::resolve_margin(masMagneticWithCoil.get_mutable_magnetic().get_mutable_coil().get_sections_description().value()[0])[0], Catch::Matchers::WithinAbs(std::max(creepageDistance, clearance) / 2, 0.00001));
        REQUIRE_THAT(OpenMagnetics::Coil::resolve_margin(masMagneticWithCoil.get_mutable_magnetic().get_mutable_coil().get_sections_description().value()[0])[1], Catch::Matchers::WithinAbs(std::max(creepageDistance, clearance) / 2, 0.00001));

        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        std::string filename = "Test_CoilAdviser_Margin_" + std::to_string(currentIndex) + ".svg";
        currentIndex++;
        outFile.append(filename);
        Painter painter(outFile);

        painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
        // painter.paint_coil_sections(masMagneticWithCoil.get_mutable_magnetic());
        painter.export_svg();
    }
    settings.reset();
}

TEST_CASE("Test_CoilAdviser_Planar", "[adviser][coil-adviser][planar]") {
    // OpenMagnetics::set_log_verbosity(2);
    settings.set_coil_maximum_layers_planar(17);
    auto gapping = OpenMagneticsTesting::get_ground_gap(0.003);
    std::vector<double> turnsRatios;
    int64_t numberStacks = 1;


    double magnetizingInductance = 10e-6;
    double temperature = 25;
    WaveformLabel waveShape = WaveformLabel::SINUSOIDAL;
    double dutyCycle = 0.5;
    double dcCurrent = 0;

    std::vector<int64_t> numberTurns = {112, 15};
    std::vector<int64_t> numberParallels = {1, 1};
    for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
        turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
    }
    double frequency = 75590;
    double peakToPeak = 13;

    auto magnetic = OpenMagneticsTesting::get_quick_magnetic("ELP 43/10/28", gapping, numberTurns, numberStacks, "3C91");
    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                          magnetizingInductance,
                                                                                          temperature,
                                                                                          waveShape,
                                                                                          peakToPeak,
                                                                                          dutyCycle,
                                                                                          dcCurrent,
                                                                                          turnsRatios);

    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(magnetic.get_mutable_core(), true);
    magnetic.get_mutable_coil().set_bobbin(bobbin);
    inputs.get_mutable_design_requirements().set_wiring_technology(MAS::WiringTechnology::PRINTED);

    OpenMagnetics::Mas masMagnetic;
    inputs.process();
    masMagnetic.set_inputs(inputs);
    masMagnetic.set_magnetic(magnetic);

    CoilAdviser coilAdviser;
    auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 100);


    REQUIRE(masMagneticsWithCoil.size() > 0);

    int currentIndex = 0;
    for (auto& masMagneticWithCoil : masMagneticsWithCoil) {
        OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        std::string filename = "Test_CoilAdviser_Planar_" + masMagneticWithCoil.get_mutable_magnetic().get_reference() + ".svg";
        currentIndex++;
        outFile.append(filename);
        Painter painter(outFile);

        painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
        // painter.paint_coil_sections(masMagneticWithCoil.get_mutable_magnetic());
        painter.export_svg();
    }
    settings.reset();
}

TEST_CASE("Test_CoilAdviser_Random_Base", "[adviser][coil-adviser][random]") {
    settings.reset();

    int count = 2;
    while (count > 0) {
        std::vector<double> turnsRatios;
        int64_t numberStacks = 1;

        std::vector<int64_t> numberTurns;
        std::vector<int64_t> numberParallels;
        std::vector<IsolationSide> isolationSides;
        int64_t numberPhysicalTurns = std::numeric_limits<int64_t>::max();
        // for (size_t windingIndex = 0; windingIndex < OpenMagnetics::TestUtils::randomSize(1, 2 + 1 - 1); ++windingIndex)
        for (size_t windingIndex = 0; windingIndex < OpenMagnetics::TestUtils::randomSize(1, 4 + 1 - 1); ++windingIndex)
        {
            int64_t numberPhysicalTurnsThisWinding = OpenMagnetics::TestUtils::randomSize(1, 300 + 1 - 1);
            int64_t numberTurnsThisWinding = OpenMagnetics::TestUtils::randomInt64(1, 100 + 1 - 1);
            numberTurns.push_back(numberTurnsThisWinding);
            isolationSides.push_back(get_isolation_side_from_index(static_cast<size_t>(OpenMagnetics::TestUtils::randomInt64(1, 10 + 1 - 1))));
            int64_t numberParallelsThisWinding = std::max(1.0, std::ceil(double(numberPhysicalTurnsThisWinding) / numberTurnsThisWinding));
            numberParallels.push_back(numberParallelsThisWinding);
            numberPhysicalTurns = std::min(numberPhysicalTurns, numberTurnsThisWinding * numberParallelsThisWinding);
        }
        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = OpenMagnetics::TestUtils::randomInt64(10000, 1000000 + 10000 - 1);
        double magnetizingInductance = double(OpenMagnetics::TestUtils::randomInt(0, 10000 - 1)) *  1e-6;
        double temperature = 25;
        WaveformLabel waveShape = WaveformLabel::SINUSOIDAL;
        double peakToPeak = OpenMagnetics::TestUtils::randomInt64(1, 30 + 1 - 1);
        double dutyCycle = double(OpenMagnetics::TestUtils::randomInt64(1, 99 + 1 - 1)) / 100;
        double dcCurrent = 0;
        if (numberTurns.size() > 1) {
            dcCurrent = OpenMagnetics::TestUtils::randomInt(0, 30 - 1);
        }
        settings.set_use_toroidal_cores(false);
        auto coreShapeNames = get_core_shape_names();
        std::string coreShapeName;
        OpenMagnetics::Magnetic magnetic;

        double gapLength = 1;
        double columnHeight = 1;
        int numberGaps = 1;
        int gappingTypeIndex = OpenMagnetics::TestUtils::randomInt(0, 4 - 1);
        GappingType gappingType = magic_enum::enum_cast<GappingType>(gappingTypeIndex).value();
        if (gappingType == GappingType::DISTRIBUTED) {
            numberGaps = OpenMagnetics::TestUtils::randomInt(0, 3 - 1);
            numberGaps *= 2;
            numberGaps += 3;
        }

        while (true) {
            coreShapeName = coreShapeNames[OpenMagnetics::TestUtils::randomSize(0, coreShapeNames.size() - 1)];
            auto shape = find_core_shape_by_name(coreShapeName);
            if (shape.get_family() == CoreShapeFamily::PQI || 
                shape.get_family() == CoreShapeFamily::UI || 
                shape.get_family() == CoreShapeFamily::EI) {
                continue;
            }

            do {
                auto shape = find_core_shape_by_name(coreShapeName);
                auto corePiece = CorePiece::factory(shape);
                gapLength = double(OpenMagnetics::TestUtils::randomInt(1.0, 10000 + static_cast<int>(1.0) - 1)) / 1000000;
                columnHeight = corePiece->get_columns()[0].get_height();
            }
            while (columnHeight < (gapLength * numberGaps));

            std::vector<CoreGap> gapping;
            switch (gappingType) {
                case GappingType::GROUND:
                    gapping = OpenMagneticsTesting::get_ground_gap(gapLength);
                    break;
                case GappingType::SPACER:
                    gapping = OpenMagneticsTesting::get_spacer_gap(gapLength);
                    break;
                case GappingType::RESIDUAL:
                    gapping = OpenMagneticsTesting::get_residual_gap();
                    break;
                case GappingType::DISTRIBUTED:
                    gapping = OpenMagneticsTesting::get_distributed_gap(gapLength, numberGaps);
                    break;
            }
            
            magnetic = OpenMagneticsTesting::get_quick_magnetic(coreShapeName, gapping, numberTurns, numberStacks, "3C91");
            break;
        }

        int turnsAlignmentIndex = OpenMagnetics::TestUtils::randomInt(0, 4 - 1);
        int windingOrientationIndex = OpenMagnetics::TestUtils::randomInt(0, 2 - 1);
        int layersOrientationIndex = OpenMagnetics::TestUtils::randomInt(0, 2 - 1);
        CoilAlignment turnsAlignment = magic_enum::enum_cast<CoilAlignment>(turnsAlignmentIndex).value();
        WindingOrientation windingOrientation = magic_enum::enum_cast<WindingOrientation>(windingOrientationIndex).value();
        WindingOrientation layersOrientation = magic_enum::enum_cast<WindingOrientation>(layersOrientationIndex).value();


        magnetic.get_mutable_coil().set_turns_alignment(turnsAlignment);
        magnetic.get_mutable_coil().set_winding_orientation(windingOrientation);
        magnetic.get_mutable_coil().set_layers_orientation(layersOrientation);

        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
        OpenMagnetics::Mas masMagnetic;
        inputs.process();
        masMagnetic.set_inputs(inputs);
        masMagnetic.set_magnetic(magnetic);

        CoilAdviser coilAdviser;
        try {
            auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);
            if (masMagneticsWithCoil.size() > 0) {
                auto masMagneticWithCoil = masMagneticsWithCoil[0];

                OpenMagneticsTesting::check_wire_standards(masMagneticWithCoil.get_mutable_magnetic().get_mutable_coil());

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
            REQUIRE(false);
            return;
        }
    }
    settings.reset();
}

TEST_CASE("Test_CoilAdviser_Random_0", "[adviser][coil-adviser][bug][smoke-test]") {
    auto gapping = OpenMagneticsTesting::get_ground_gap(0.003);
    std::vector<double> turnsRatios;
    int64_t numberStacks = 1;


    double magnetizingInductance = 10e-6;
    double temperature = 25;
    WaveformLabel waveShape = WaveformLabel::SINUSOIDAL;
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
    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                          magnetizingInductance,
                                                                                          temperature,
                                                                                          waveShape,
                                                                                          peakToPeak,
                                                                                          dutyCycle,
                                                                                          dcCurrent,
                                                                                          turnsRatios);

    OpenMagnetics::Mas masMagnetic;
    inputs.process();
    masMagnetic.set_inputs(inputs);
    masMagnetic.set_magnetic(magnetic);

    CoilAdviser coilAdviser;
    auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

    if (masMagneticsWithCoil.size() > 0) {
        auto masMagneticWithCoil = masMagneticsWithCoil[0];
        OpenMagneticsTesting::check_wire_standards(masMagneticWithCoil.get_mutable_magnetic().get_mutable_coil());
        OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
    }
}

TEST_CASE("Test_CoilAdviser_Random_1", "[adviser][coil-adviser][bug][smoke-test]") {
    std::vector<double> turnsRatios;
    int64_t numberStacks = 1;


    double magnetizingInductance = 10e-6;
    double temperature = 25;
    WaveformLabel waveShape = WaveformLabel::SINUSOIDAL;
    double dutyCycle = 0.5;
    double dcCurrent = 0;

    std::vector<int64_t> numberTurns = {16, 34};
    std::vector<int64_t> numberParallels = {1, 1};
    for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
        turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
    }
    double frequency = 811022;
    double peakToPeak = 1;

    OpenMagnetics::Magnetic magnetic;
    std::string coreShapeName  = "P 7.35X3.6";
    double gapLength = 1;
    double columnHeight = 1;
    do {
        auto shape = find_core_shape_by_name(coreShapeName);
        auto corePiece = CorePiece::factory(shape);
        gapLength = double(OpenMagnetics::TestUtils::randomInt(1.0, 10000 + static_cast<int>(1.0) - 1)) / 1000000;
        columnHeight = corePiece->get_columns()[0].get_height();
    }
    while (columnHeight < gapLength);
    auto gapping = OpenMagneticsTesting::get_ground_gap(gapLength);
    magnetic = OpenMagneticsTesting::get_quick_magnetic(coreShapeName, gapping, numberTurns, numberStacks, "3C91");
    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                          magnetizingInductance,
                                                                                          temperature,
                                                                                          waveShape,
                                                                                          peakToPeak,
                                                                                          dutyCycle,
                                                                                          dcCurrent,
                                                                                          turnsRatios);

    OpenMagnetics::Mas masMagnetic;
    inputs.process();
    masMagnetic.set_inputs(inputs);
    masMagnetic.set_magnetic(magnetic);

    CoilAdviser coilAdviser;
    auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

    if (masMagneticsWithCoil.size() > 0) {
        auto masMagneticWithCoil = masMagneticsWithCoil[0];
        OpenMagneticsTesting::check_wire_standards(masMagneticWithCoil.get_mutable_magnetic().get_mutable_coil());
        OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        std::string filename = "Test_CoilAdviser" + std::to_string(OpenMagnetics::TestUtils::randomInt(0, RAND_MAX)) + ".svg";
        outFile.append(filename);
        // Painter painter(outFile);

        // painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
        // painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
        // painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
        // painter.export_svg();
    }
}

TEST_CASE("Test_CoilAdviser_Random_2", "[adviser][coil-adviser][bug]") {

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
    WaveformLabel waveShape = WaveformLabel::SINUSOIDAL;
    double peakToPeak = 10;
    double dutyCycle = 0.5;
    double dcCurrent = 0;
    settings.set_use_toroidal_cores(false);
    auto coreShapeNames = get_core_shape_names();
    std::string coreShapeName;
    OpenMagnetics::Magnetic magnetic;

    auto gapping = OpenMagneticsTesting::get_residual_gap();
    magnetic = OpenMagneticsTesting::get_quick_magnetic("ETD 19", gapping, numberTurns, numberStacks, "3C91");


    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                          magnetizingInductance,
                                                                                          temperature,
                                                                                          waveShape,
                                                                                          peakToPeak,
                                                                                          dutyCycle,
                                                                                          dcCurrent,
                                                                                          turnsRatios);

    OpenMagnetics::Mas masMagnetic;
    inputs.process();
    masMagnetic.set_inputs(inputs);
    masMagnetic.set_magnetic(magnetic);

    CoilAdviser coilAdviser;
    auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

    if (masMagneticsWithCoil.size() > 0) {
        auto masMagneticWithCoil = masMagneticsWithCoil[0];
        OpenMagneticsTesting::check_wire_standards(masMagneticWithCoil.get_mutable_magnetic().get_mutable_coil());
        OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        std::string filename = "Test_CoilAdviser" + std::to_string(OpenMagnetics::TestUtils::randomInt(0, RAND_MAX)) + ".svg";
        outFile.append(filename);
        // Painter painter(outFile);

        // painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
        // painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
        // painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
        // painter.export_svg();
    }
}

TEST_CASE("Test_CoilAdviser_Random_3", "[adviser][coil-adviser][bug][smoke-test]") {

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
    WaveformLabel waveShape = WaveformLabel::SINUSOIDAL;
    double peakToPeak = 18;
    double dutyCycle = 0.88;
    double dcCurrent = 6;
    settings.set_use_toroidal_cores(false);
    auto coreShapeNames = get_core_shape_names();
    std::string coreShapeName;
    OpenMagnetics::Magnetic magnetic;

    CoilAlignment turnsAlignment = magic_enum::enum_cast<CoilAlignment>(3).value();
    WindingOrientation windingOrientation = magic_enum::enum_cast<WindingOrientation>(0).value();
    WindingOrientation layersOrientation = magic_enum::enum_cast<WindingOrientation>(1).value();
    auto gapping = OpenMagneticsTesting::get_ground_gap(0.00756);
    magnetic = OpenMagneticsTesting::get_quick_magnetic("U 81/39/20", gapping, numberTurns, numberStacks, "3C91");

    magnetic.get_mutable_coil().set_turns_alignment(turnsAlignment);
    magnetic.get_mutable_coil().set_winding_orientation(windingOrientation);
    magnetic.get_mutable_coil().set_layers_orientation(layersOrientation);

    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                          magnetizingInductance,
                                                                                          temperature,
                                                                                          waveShape,
                                                                                          peakToPeak,
                                                                                          dutyCycle,
                                                                                          dcCurrent,
                                                                                          turnsRatios);

    OpenMagnetics::Mas masMagnetic;
    inputs.process();
    masMagnetic.set_inputs(inputs);
    masMagnetic.set_magnetic(magnetic);

    CoilAdviser coilAdviser;
    auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

    if (masMagneticsWithCoil.size() > 0) {
        auto masMagneticWithCoil = masMagneticsWithCoil[0];
        OpenMagneticsTesting::check_wire_standards(masMagneticWithCoil.get_mutable_magnetic().get_mutable_coil());
        OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        std::string filename = "Test_CoilAdviser" + std::to_string(OpenMagnetics::TestUtils::randomInt(0, RAND_MAX)) + ".svg";
        outFile.append(filename);
        Painter painter(outFile);

        painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
        painter.export_svg();
    }
}

TEST_CASE("Test_CoilAdviser_Random_4", "[adviser][coil-adviser][bug][smoke-test]") {
    settings.set_coil_wind_even_if_not_fit(false);

    std::vector<double> turnsRatios;
    int64_t numberStacks = 1;

    std::vector<int64_t> numberTurns = {28};
    std::vector<int64_t> numberParallels = {76};
    std::vector<IsolationSide> isolationSides = {IsolationSide::OCTONARY, IsolationSide::QUINARY};

    double frequency = 837961;
    double magnetizingInductance = 0.007191;
    double temperature = 25;
    WaveformLabel waveShape = WaveformLabel::SINUSOIDAL;
    double peakToPeak = 3;
    double dutyCycle = 0.18;
    double dcCurrent = 28;
    settings.set_use_toroidal_cores(false);
    auto coreShapeNames = get_core_shape_names();
    std::string coreShapeName;
    OpenMagnetics::Magnetic magnetic;

    CoilAlignment turnsAlignment = magic_enum::enum_cast<CoilAlignment>(2).value();
    WindingOrientation windingOrientation = magic_enum::enum_cast<WindingOrientation>(0).value();
    WindingOrientation layersOrientation = magic_enum::enum_cast<WindingOrientation>(0).value();
    auto gapping = OpenMagneticsTesting::get_distributed_gap(0.004025, 3);
    magnetic = OpenMagneticsTesting::get_quick_magnetic("UR 39/35/15", gapping, numberTurns, numberStacks, "3C91");

    magnetic.get_mutable_coil().set_turns_alignment(turnsAlignment);
    magnetic.get_mutable_coil().set_winding_orientation(windingOrientation);
    magnetic.get_mutable_coil().set_layers_orientation(layersOrientation);

    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                          magnetizingInductance,
                                                                                          temperature,
                                                                                          waveShape,
                                                                                          peakToPeak,
                                                                                          dutyCycle,
                                                                                          dcCurrent,
                                                                                          turnsRatios);

    inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
    OpenMagnetics::Mas masMagnetic;
    inputs.process();
    masMagnetic.set_inputs(inputs);
    masMagnetic.set_magnetic(magnetic);

    CoilAdviser coilAdviser;
    auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

    if (masMagneticsWithCoil.size() > 0) {
        auto masMagneticWithCoil = masMagneticsWithCoil[0];
        OpenMagneticsTesting::check_wire_standards(masMagneticWithCoil.get_mutable_magnetic().get_mutable_coil());
        OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        std::string filename = "Test_CoilAdviser" + std::to_string(OpenMagnetics::TestUtils::randomInt(0, RAND_MAX)) + ".svg";
        outFile.append(filename);
        Painter painter(outFile);

        painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
        painter.export_svg();
    }
    settings.reset();
}

TEST_CASE("Test_CoilAdviser_Random_5", "[adviser][coil-adviser][bug]") {

    int64_t numberStacks = 1;

    std::vector<int64_t> numberTurns = {15, 36, 87, 60};
    std::vector<int64_t> numberParallels = {1, 1, 1, 1};
    std::vector<double> turnsRatios = {15. / 36, 15. / 87, 15. / 60};
    std::vector<IsolationSide> isolationSides = {IsolationSide::DENARY,
                                                                IsolationSide::NONARY,
                                                                IsolationSide::QUATERNARY,
                                                                IsolationSide::OCTONARY};

    double frequency = 592535;
    double magnetizingInductance = 0.002575;
    double temperature = 23;
    WaveformLabel waveShape = WaveformLabel::SINUSOIDAL;
    double peakToPeak = 3;
    double dutyCycle = 0.61;
    double dcCurrent = 5;
    settings.set_use_toroidal_cores(false);
    auto coreShapeNames = get_core_shape_names();
    std::string coreShapeName;
    OpenMagnetics::Magnetic magnetic;

    CoilAlignment turnsAlignment = magic_enum::enum_cast<CoilAlignment>(3).value();
    WindingOrientation windingOrientation = magic_enum::enum_cast<WindingOrientation>(0).value();
    WindingOrientation layersOrientation = magic_enum::enum_cast<WindingOrientation>(1).value();
    auto gapping = OpenMagneticsTesting::get_distributed_gap(0.003175, 1);
    magnetic = OpenMagneticsTesting::get_quick_magnetic("E 114/46/26", gapping, numberTurns, numberStacks, "3C91");

    magnetic.get_mutable_coil().set_turns_alignment(turnsAlignment);
    magnetic.get_mutable_coil().set_winding_orientation(windingOrientation);
    magnetic.get_mutable_coil().set_layers_orientation(layersOrientation);

    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                          magnetizingInductance,
                                                                                          temperature,
                                                                                          waveShape,
                                                                                          peakToPeak,
                                                                                          dutyCycle,
                                                                                          dcCurrent,
                                                                                          turnsRatios);

    inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
    OpenMagnetics::Mas masMagnetic;
    inputs.process();
    masMagnetic.set_inputs(inputs);
    masMagnetic.set_magnetic(magnetic);

    CoilAdviser coilAdviser;
    auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

    if (masMagneticsWithCoil.size() > 0) {
        auto masMagneticWithCoil = masMagneticsWithCoil[0];
        OpenMagneticsTesting::check_wire_standards(masMagneticWithCoil.get_mutable_magnetic().get_mutable_coil());
        OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        std::string filename = "Test_CoilAdviser" + std::to_string(OpenMagnetics::TestUtils::randomInt(0, RAND_MAX)) + ".svg";
        outFile.append(filename);
        Painter painter(outFile);

        painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
        painter.export_svg();
    }
}

TEST_CASE("Test_CoilAdviser_Random_6", "[adviser][coil-adviser][bug][smoke-test]") {

    int64_t numberStacks = 1;

    std::vector<int64_t> numberTurns = {11, 82};
    std::vector<int64_t> numberParallels = {3, 3};
    std::vector<double> turnsRatios = {11. / 82};
    std::vector<IsolationSide> isolationSides = {IsolationSide::DENARY,
                                                                IsolationSide::NONARY};

    double frequency = 617645;
    double magnetizingInductance = 0.009088;
    double temperature = 23;
    WaveformLabel waveShape = WaveformLabel::SINUSOIDAL;
    double peakToPeak = 26;
    double dutyCycle = 0.38;
    double dcCurrent = 16;
    settings.set_use_toroidal_cores(false);
    auto coreShapeNames = get_core_shape_names();
    std::string coreShapeName;
    OpenMagnetics::Magnetic magnetic;

    CoilAlignment turnsAlignment = magic_enum::enum_cast<CoilAlignment>(1).value();
    WindingOrientation windingOrientation = magic_enum::enum_cast<WindingOrientation>(0).value();
    WindingOrientation layersOrientation = magic_enum::enum_cast<WindingOrientation>(1).value();
    auto gapping = OpenMagneticsTesting::get_ground_gap(0.001227);
    magnetic = OpenMagneticsTesting::get_quick_magnetic("E 72/28/19", gapping, numberTurns, numberStacks, "3C91");

    magnetic.get_mutable_coil().set_turns_alignment(turnsAlignment);
    magnetic.get_mutable_coil().set_winding_orientation(windingOrientation);
    magnetic.get_mutable_coil().set_layers_orientation(layersOrientation);

    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                          magnetizingInductance,
                                                                                          temperature,
                                                                                          waveShape,
                                                                                          peakToPeak,
                                                                                          dutyCycle,
                                                                                          dcCurrent,
                                                                                          turnsRatios);

    inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
    OpenMagnetics::Mas masMagnetic;
    inputs.process();
    masMagnetic.set_inputs(inputs);
    masMagnetic.set_magnetic(magnetic);

    CoilAdviser coilAdviser;
    auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

    if (masMagneticsWithCoil.size() > 0) {
        auto masMagneticWithCoil = masMagneticsWithCoil[0];
        OpenMagneticsTesting::check_wire_standards(masMagneticWithCoil.get_mutable_magnetic().get_mutable_coil());
        OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        std::string filename = "Test_CoilAdviser" + std::to_string(OpenMagnetics::TestUtils::randomInt(0, RAND_MAX)) + ".svg";
        outFile.append(filename);
        Painter painter(outFile);

        painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
        painter.export_svg();
    }
}

TEST_CASE("Test_CoilAdviser_Random_7", "[adviser][coil-adviser][bug]") {

    int64_t numberStacks = 1;

    std::vector<int64_t> numberTurns = {60, 18, 10};
    std::vector<int64_t> numberParallels = {1, 3, 3};
    std::vector<double> turnsRatios = {60. / 18, 60. / 10};
    std::vector<IsolationSide> isolationSides = {IsolationSide::SECONDARY,
                                                                IsolationSide::NONARY,
                                                                IsolationSide::TERTIARY};

    double frequency = 95989;
    double magnetizingInductance = 0.009266;
    double temperature = 23;
    WaveformLabel waveShape = WaveformLabel::SINUSOIDAL;
    double peakToPeak = 9;
    double dutyCycle = 0.57;
    double dcCurrent = 1;
    settings.set_use_toroidal_cores(false);
    auto coreShapeNames = get_core_shape_names();
    std::string coreShapeName;
    OpenMagnetics::Magnetic magnetic;

    CoilAlignment turnsAlignment = magic_enum::enum_cast<CoilAlignment>(1).value();
    WindingOrientation windingOrientation = magic_enum::enum_cast<WindingOrientation>(1).value();
    WindingOrientation layersOrientation = magic_enum::enum_cast<WindingOrientation>(0).value();
    auto gapping = OpenMagneticsTesting::get_ground_gap(0.007023);
    magnetic = OpenMagneticsTesting::get_quick_magnetic("E 35", gapping, numberTurns, numberStacks, "3C91");

    magnetic.get_mutable_coil().set_turns_alignment(turnsAlignment);
    magnetic.get_mutable_coil().set_winding_orientation(windingOrientation);
    magnetic.get_mutable_coil().set_layers_orientation(layersOrientation);

    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                          magnetizingInductance,
                                                                                          temperature,
                                                                                          waveShape,
                                                                                          peakToPeak,
                                                                                          dutyCycle,
                                                                                          dcCurrent,
                                                                                          turnsRatios);

    inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
    OpenMagnetics::Mas masMagnetic;
    inputs.process();
    masMagnetic.set_inputs(inputs);
    masMagnetic.set_magnetic(magnetic);

    CoilAdviser coilAdviser;
    auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

    if (masMagneticsWithCoil.size() > 0) {
        auto masMagneticWithCoil = masMagneticsWithCoil[0];
        OpenMagneticsTesting::check_wire_standards(masMagneticWithCoil.get_mutable_magnetic().get_mutable_coil());
        OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        std::string filename = "Test_CoilAdviser" + std::to_string(OpenMagnetics::TestUtils::randomInt(0, RAND_MAX)) + ".svg";
        outFile.append(filename);
        Painter painter(outFile);

        painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
        painter.export_svg();
    }
}

TEST_CASE("Test_CoilAdviser_Random_8", "[adviser][coil-adviser][bug]") {

    int64_t numberStacks = 1;

    std::vector<int64_t> numberTurns = {36, 3, 8};
    std::vector<int64_t> numberParallels = {2, 1, 2};
    std::vector<double> turnsRatios = {36. / 3, 36. / 8};
    std::vector<IsolationSide> isolationSides = {IsolationSide::UNDENARY,
                                                                IsolationSide::OCTONARY,
                                                                IsolationSide::DENARY};

    double frequency = 632226;
    double magnetizingInductance = 0.001529;
    double temperature = 23;
    WaveformLabel waveShape = WaveformLabel::SINUSOIDAL;
    double peakToPeak = 3;
    double dutyCycle = 0.68;
    double dcCurrent = 7;
    settings.set_use_toroidal_cores(false);
    auto coreShapeNames = get_core_shape_names();
    std::string coreShapeName;
    OpenMagnetics::Magnetic magnetic;

    CoilAlignment turnsAlignment = magic_enum::enum_cast<CoilAlignment>(0).value();
    WindingOrientation windingOrientation = magic_enum::enum_cast<WindingOrientation>(0).value();
    WindingOrientation layersOrientation = magic_enum::enum_cast<WindingOrientation>(0).value();
    auto gapping = OpenMagneticsTesting::get_ground_gap(0.009828);
    magnetic = OpenMagneticsTesting::get_quick_magnetic("E 80/38/20", gapping, numberTurns, numberStacks, "3C91");

    magnetic.get_mutable_coil().set_turns_alignment(turnsAlignment);
    magnetic.get_mutable_coil().set_winding_orientation(windingOrientation);
    magnetic.get_mutable_coil().set_layers_orientation(layersOrientation);

    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                          magnetizingInductance,
                                                                                          temperature,
                                                                                          waveShape,
                                                                                          peakToPeak,
                                                                                          dutyCycle,
                                                                                          dcCurrent,
                                                                                          turnsRatios);

    inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
    OpenMagnetics::Mas masMagnetic;
    inputs.process();
    masMagnetic.set_inputs(inputs);
    masMagnetic.set_magnetic(magnetic);

    CoilAdviser coilAdviser;
    auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

    if (masMagneticsWithCoil.size() > 0) {
        auto masMagneticWithCoil = masMagneticsWithCoil[0];
        OpenMagneticsTesting::check_wire_standards(masMagneticWithCoil.get_mutable_magnetic().get_mutable_coil());
        OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        std::string filename = "Test_CoilAdviser" + std::to_string(OpenMagnetics::TestUtils::randomInt(0, RAND_MAX)) + ".svg";
        outFile.append(filename);
        Painter painter(outFile);

        painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
        painter.export_svg();
    }
}

TEST_CASE("Test_CoilAdviser_Random_9", "[adviser][coil-adviser][bug]") {

    int64_t numberStacks = 1;

    std::vector<int64_t> numberTurns = {36, 55, 96};
    std::vector<int64_t> numberParallels = {2, 1, 2};
    std::vector<double> turnsRatios = {36. / 55, 36. / 96};
    std::vector<IsolationSide> isolationSides = {IsolationSide::SENARY,
                                                                IsolationSide::NONARY,
                                                                IsolationSide::NONARY};

    double frequency = 632226;
    double magnetizingInductance = 0.001529;
    double temperature = 23;
    WaveformLabel waveShape = WaveformLabel::SINUSOIDAL;
    double peakToPeak = 3;
    double dutyCycle = 0.68;
    double dcCurrent = 7;
    settings.set_use_toroidal_cores(false);
    auto coreShapeNames = get_core_shape_names();
    std::string coreShapeName;
    OpenMagnetics::Magnetic magnetic;

    CoilAlignment turnsAlignment = magic_enum::enum_cast<CoilAlignment>(2).value();
    WindingOrientation windingOrientation = magic_enum::enum_cast<WindingOrientation>(0).value();
    WindingOrientation layersOrientation = magic_enum::enum_cast<WindingOrientation>(0).value();
    auto gapping = OpenMagneticsTesting::get_ground_gap(0.009828);
    magnetic = OpenMagneticsTesting::get_quick_magnetic("U 79/129/31", gapping, numberTurns, numberStacks, "3C91");

    magnetic.get_mutable_coil().set_turns_alignment(turnsAlignment);
    magnetic.get_mutable_coil().set_winding_orientation(windingOrientation);
    magnetic.get_mutable_coil().set_layers_orientation(layersOrientation);

    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                          magnetizingInductance,
                                                                                          temperature,
                                                                                          waveShape,
                                                                                          peakToPeak,
                                                                                          dutyCycle,
                                                                                          dcCurrent,
                                                                                          turnsRatios);

    inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
    OpenMagnetics::Mas masMagnetic;
    inputs.process();
    masMagnetic.set_inputs(inputs);
    masMagnetic.set_magnetic(magnetic);

    CoilAdviser coilAdviser;
    auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

    if (masMagneticsWithCoil.size() > 0) {
        auto masMagneticWithCoil = masMagneticsWithCoil[0];
        OpenMagneticsTesting::check_wire_standards(masMagneticWithCoil.get_mutable_magnetic().get_mutable_coil());
        OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        std::string filename = "Test_CoilAdviser" + std::to_string(OpenMagnetics::TestUtils::randomInt(0, RAND_MAX)) + ".svg";
        outFile.append(filename);
        Painter painter(outFile);

        painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
        painter.export_svg();
    }
}

TEST_CASE("Test_CoilAdviser_Random_10", "[adviser][coil-adviser][bug]") {

    int64_t numberStacks = 1;

    std::vector<int64_t> numberTurns = {49, 80, 78, 1};
    std::vector<int64_t> numberParallels = {3, 1, 1, 2};
    std::vector<double> turnsRatios = {49. / 80, 49. / 78, 49};
    std::vector<IsolationSide> isolationSides = {IsolationSide::TERTIARY,
                                                                IsolationSide::SENARY,
                                                                IsolationSide::DENARY,
                                                                IsolationSide::NONARY};

    double frequency = 660462;
    double magnetizingInductance = 0.006606;
    double temperature = 23;
    WaveformLabel waveShape = WaveformLabel::SINUSOIDAL;
    double peakToPeak = 24;
    double dutyCycle = 0.28;
    double dcCurrent = 2;
    settings.set_use_toroidal_cores(false);
    auto coreShapeNames = get_core_shape_names();
    std::string coreShapeName;
    OpenMagnetics::Magnetic magnetic;

    CoilAlignment turnsAlignment = magic_enum::enum_cast<CoilAlignment>(2).value();
    WindingOrientation windingOrientation = magic_enum::enum_cast<WindingOrientation>(1).value();
    WindingOrientation layersOrientation = magic_enum::enum_cast<WindingOrientation>(1).value();
    auto gapping = OpenMagneticsTesting::get_distributed_gap(0.002048, 3);
    magnetic = OpenMagneticsTesting::get_quick_magnetic("E 50/15", gapping, numberTurns, numberStacks, "3C91");

    magnetic.get_mutable_coil().set_turns_alignment(turnsAlignment);
    magnetic.get_mutable_coil().set_winding_orientation(windingOrientation);
    magnetic.get_mutable_coil().set_layers_orientation(layersOrientation);

    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                          magnetizingInductance,
                                                                                          temperature,
                                                                                          waveShape,
                                                                                          peakToPeak,
                                                                                          dutyCycle,
                                                                                          dcCurrent,
                                                                                          turnsRatios);

    inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
    OpenMagnetics::Mas masMagnetic;
    inputs.process();
    masMagnetic.set_inputs(inputs);
    masMagnetic.set_magnetic(magnetic);

    CoilAdviser coilAdviser;
    auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

    if (masMagneticsWithCoil.size() > 0) {
        auto masMagneticWithCoil = masMagneticsWithCoil[0];
        // OpenMagneticsTesting::check_wire_standards(masMagneticWithCoil.get_mutable_magnetic().get_mutable_coil());
        // OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        std::string filename = "Test_CoilAdviser" + std::to_string(OpenMagnetics::TestUtils::randomInt(0, RAND_MAX)) + ".svg";
        outFile.append(filename);
        Painter painter(outFile);

        painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
        painter.export_svg();
    }
}

TEST_CASE("Test_CoilAdviser_Random_11", "[adviser][coil-adviser][bug][smoke-test]") {

    int64_t numberStacks = 1;

    std::vector<int64_t> numberTurns = {72};
    std::vector<int64_t> numberParallels = {4};
    std::vector<double> turnsRatios = {};
    std::vector<IsolationSide> isolationSides = {IsolationSide::NONARY};

    double frequency = 821021;
    double magnetizingInductance = 8.6e-05;
    double temperature = 23;
    WaveformLabel waveShape = WaveformLabel::SINUSOIDAL;
    double peakToPeak = 19;
    double dutyCycle = 0.14;
    double dcCurrent = 0;
    settings.set_use_toroidal_cores(false);
    auto coreShapeNames = get_core_shape_names();
    std::string coreShapeName;
    OpenMagnetics::Magnetic magnetic;

    CoilAlignment turnsAlignment = magic_enum::enum_cast<CoilAlignment>(3).value();
    WindingOrientation windingOrientation = magic_enum::enum_cast<WindingOrientation>(0).value();
    WindingOrientation layersOrientation = magic_enum::enum_cast<WindingOrientation>(1).value();
    auto gapping = OpenMagneticsTesting::get_distributed_gap(0.001424, 1);
    magnetic = OpenMagneticsTesting::get_quick_magnetic("U 15/11/6", gapping, numberTurns, numberStacks, "3C91");

    magnetic.get_mutable_coil().set_turns_alignment(turnsAlignment);
    magnetic.get_mutable_coil().set_winding_orientation(windingOrientation);
    magnetic.get_mutable_coil().set_layers_orientation(layersOrientation);

    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                          magnetizingInductance,
                                                                                          temperature,
                                                                                          waveShape,
                                                                                          peakToPeak,
                                                                                          dutyCycle,
                                                                                          dcCurrent,
                                                                                          turnsRatios);

    inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
    OpenMagnetics::Mas masMagnetic;
    inputs.process();
    masMagnetic.set_inputs(inputs);
    masMagnetic.set_magnetic(magnetic);

    CoilAdviser coilAdviser;
    auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

    if (masMagneticsWithCoil.size() > 0) {
        auto masMagneticWithCoil = masMagneticsWithCoil[0];
        // OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        std::string filename = "Test_CoilAdviser" + std::to_string(OpenMagnetics::TestUtils::randomInt(0, RAND_MAX)) + ".svg";
        outFile.append(filename);
        Painter painter(outFile);

        painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
        painter.export_svg();
    }
}

TEST_CASE("Test_CoilAdviser_Random_12", "[adviser][coil-adviser][bug]") {

    int64_t numberStacks = 1;

    std::vector<int64_t> numberTurns = {53, 100, 80, 98};
    std::vector<int64_t> numberParallels = {5, 2, 3, 2};
    std::vector<double> turnsRatios = {53.0 / 100, 53.0 / 80, 53.0 / 98};
    std::vector<IsolationSide> isolationSides = {IsolationSide::OCTONARY,
                                                                IsolationSide::SENARY,
                                                                IsolationSide::QUINARY,
                                                                IsolationSide::SENARY};


    double frequency = 460425;
    double magnetizingInductance = 0.005275;
    double temperature = 23;
    WaveformLabel waveShape = WaveformLabel::SINUSOIDAL;
    double peakToPeak = 28;
    double dutyCycle = 0.73;
    double dcCurrent = 5;
    settings.set_use_toroidal_cores(false);
    auto coreShapeNames = get_core_shape_names();
    std::string coreShapeName;
    OpenMagnetics::Magnetic magnetic;

    CoilAlignment turnsAlignment = magic_enum::enum_cast<CoilAlignment>(0).value();
    WindingOrientation windingOrientation = magic_enum::enum_cast<WindingOrientation>(1).value();
    WindingOrientation layersOrientation = magic_enum::enum_cast<WindingOrientation>(0).value();

    auto gapping = OpenMagneticsTesting::get_spacer_gap(0.006456);
    magnetic = OpenMagneticsTesting::get_quick_magnetic("U 81/39/20", gapping, numberTurns, numberStacks, "3C91");

    magnetic.get_mutable_coil().set_turns_alignment(turnsAlignment);
    magnetic.get_mutable_coil().set_winding_orientation(windingOrientation);
    magnetic.get_mutable_coil().set_layers_orientation(layersOrientation);

    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                          magnetizingInductance,
                                                                                          temperature,
                                                                                          waveShape,
                                                                                          peakToPeak,
                                                                                          dutyCycle,
                                                                                          dcCurrent,
                                                                                          turnsRatios);

    inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
    OpenMagnetics::Mas masMagnetic;
    inputs.process();
    masMagnetic.set_inputs(inputs);
    masMagnetic.set_magnetic(magnetic);

    CoilAdviser coilAdviser;
    auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

    if (masMagneticsWithCoil.size() > 0) {
        auto masMagneticWithCoil = masMagneticsWithCoil[0];
        OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        std::string filename = "Test_CoilAdviser_Random_12.svg";
        outFile.append(filename);
        Painter painter(outFile);

        painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
        painter.export_svg();
    }
}

TEST_CASE("Test_CoilAdviser_Random_13", "[adviser][coil-adviser][bug][smoke-test]") {

    int64_t numberStacks = 1;

    std::vector<int64_t> numberTurns = {63};
    std::vector<int64_t> numberParallels = {3};
    std::vector<double> turnsRatios = {};
    std::vector<IsolationSide> isolationSides = {IsolationSide::NONARY};


    double frequency = 334200;
    double peakToPeak = 13;
    double magnetizingInductance = 0.008387;
    double dutyCycle = 0.33;
    double dcCurrent = 0;
    double temperature = 23;
    WaveformLabel waveShape = WaveformLabel::SINUSOIDAL;
    settings.set_use_toroidal_cores(false);
    auto coreShapeNames = get_core_shape_names();
    std::string coreShapeName;
    OpenMagnetics::Magnetic magnetic;

    CoilAlignment turnsAlignment = magic_enum::enum_cast<CoilAlignment>(1).value();
    WindingOrientation windingOrientation = magic_enum::enum_cast<WindingOrientation>(1).value();
    WindingOrientation layersOrientation = magic_enum::enum_cast<WindingOrientation>(0).value();

    auto gapping = OpenMagneticsTesting::get_spacer_gap(0.002261);
    magnetic = OpenMagneticsTesting::get_quick_magnetic("EQ 32/22/18", gapping, numberTurns, numberStacks, "3C91");

    magnetic.get_mutable_coil().set_turns_alignment(turnsAlignment);
    magnetic.get_mutable_coil().set_winding_orientation(windingOrientation);
    magnetic.get_mutable_coil().set_layers_orientation(layersOrientation);

    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                          magnetizingInductance,
                                                                                          temperature,
                                                                                          waveShape,
                                                                                          peakToPeak,
                                                                                          dutyCycle,
                                                                                          dcCurrent,
                                                                                          turnsRatios);

    inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
    OpenMagnetics::Mas masMagnetic;
    inputs.process();
    masMagnetic.set_inputs(inputs);
    masMagnetic.set_magnetic(magnetic);

    CoilAdviser coilAdviser;
    auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

    if (masMagneticsWithCoil.size() > 0) {
        auto masMagneticWithCoil = masMagneticsWithCoil[0];
        OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        std::string filename = "Test_CoilAdviser_Random_13.svg";
        outFile.append(filename);
        Painter painter(outFile);

        painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
        painter.export_svg();
    }
}

TEST_CASE("Test_CoilAdviser_Random_14", "[adviser][coil-adviser][bug][smoke-test]") {

    int64_t numberStacks = 1;

    std::vector<int64_t> numberTurns = {73, 32, 1};
    std::vector<int64_t> numberParallels = {1, 10, 226};
    std::vector<double> turnsRatios = {32.0 / 73, 1.0 / 73};
    std::vector<IsolationSide> isolationSides = {IsolationSide::OCTONARY,
                                                                IsolationSide::DENARY,
                                                                IsolationSide::TERTIARY};


    double frequency = 53298;
    double peakToPeak = 13;
    double magnetizingInductance = 0.005401;
    double dutyCycle = 0.93;
    double dcCurrent = 8;
    double temperature = 23;
    WaveformLabel waveShape = WaveformLabel::SINUSOIDAL;
    settings.set_use_toroidal_cores(false);
    auto coreShapeNames = get_core_shape_names();
    std::string coreShapeName;
    OpenMagnetics::Magnetic magnetic;

    CoilAlignment turnsAlignment = magic_enum::enum_cast<CoilAlignment>(0).value();
    WindingOrientation windingOrientation = magic_enum::enum_cast<WindingOrientation>(1).value();
    WindingOrientation layersOrientation = magic_enum::enum_cast<WindingOrientation>(1).value();

    auto gapping = OpenMagneticsTesting::get_ground_gap(0.000503);
    magnetic = OpenMagneticsTesting::get_quick_magnetic("UR 35/27.5/13", gapping, numberTurns, numberStacks, "3C91");

    magnetic.get_mutable_coil().set_turns_alignment(turnsAlignment);
    magnetic.get_mutable_coil().set_winding_orientation(windingOrientation);
    magnetic.get_mutable_coil().set_layers_orientation(layersOrientation);

    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                          magnetizingInductance,
                                                                                          temperature,
                                                                                          waveShape,
                                                                                          peakToPeak,
                                                                                          dutyCycle,
                                                                                          dcCurrent,
                                                                                          turnsRatios);

    inputs.get_mutable_design_requirements().set_isolation_sides(isolationSides);
    OpenMagnetics::Mas masMagnetic;
    inputs.process();
    masMagnetic.set_inputs(inputs);
    masMagnetic.set_magnetic(magnetic);

    CoilAdviser coilAdviser;
    auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

    if (masMagneticsWithCoil.size() > 0) {
        auto masMagneticWithCoil = masMagneticsWithCoil[0];
        OpenMagneticsTesting::check_wire_standards(masMagneticWithCoil.get_mutable_magnetic().get_mutable_coil());
        OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        std::string filename = "Test_CoilAdviser_Random_14.svg";
        outFile.append(filename);
        Painter painter(outFile);

        painter.paint_core(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_bobbin(masMagneticWithCoil.get_mutable_magnetic());
        painter.paint_coil_turns(masMagneticWithCoil.get_mutable_magnetic());
        painter.export_svg();
    }
}

TEST_CASE("Test_CoilAdviser_Random_15", "[adviser][coil-adviser][bug][smoke-test]") {
    OpenMagnetics::Mas mas(json::parse(R"({"inputs":{"designRequirements":{"magnetizingInductance":{"nominal":0.0001},"minimumImpedance":[{"frequency":100000,"impedance":{"magnitude":1000}}],"name":"My Design Requirements","turnsRatios":[{"nominal":1}]},"operatingPoints":[{"name":"Operating Point No. 1","conditions":{"ambientTemperature":42},"excitationsPerWinding":[{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"data":[-5,5,-5],"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"Triangular","acEffectiveFrequency":110746.40291779551,"effectiveFrequency":110746.40291779551,"peak":5,"rms":2.8874560332150576,"thd":0.12151487440704967},"harmonics":{"amplitudes":[1.1608769501236793e-14,4.05366124583194,1.787369544444173e-15,0.4511310569983995,9.749053004706756e-16,0.16293015292554872,4.036157626725542e-16,0.08352979924600704,3.4998295008010614e-16,0.0508569581336163,3.1489164048780735e-16,0.034320410449418075,3.142469873118059e-16,0.024811988673843106,2.3653352035940994e-16,0.018849001010678823,2.9306524147249266e-16,0.014866633059596499,1.796485796132569e-16,0.012077180559557785,1.6247782523152451e-16,0.010049063750920609,1.5324769149805092e-16,0.008529750975091871,1.0558579733068502e-16,0.007363501410705499,7.513269775674661e-17,0.006450045785294609,5.871414177162291e-17,0.005722473794997712,9.294731722001391e-17,0.005134763398167541,1.194820309200107e-16,0.004654430423785411,8.2422739080512e-17,0.004258029771397032,9.5067306351894e-17,0.0039283108282380024,1.7540347128474968e-16,0.0036523670873925395,9.623794010508822e-17,0.0034204021424253787,1.4083470894369491e-16,0.0032248884817922927,1.4749333016985644e-16,0.0030599828465501895,1.0448590642474364e-16,0.002921112944200383,7.575487373767413e-18,0.002804680975178716,7.419510610361002e-17,0.0027078483284668376,3.924741709073613e-17,0.0026283777262804953,2.2684279102637236e-17,0.0025645167846443107,8.997077625295079e-17,0.0025149120164513483,7.131074184849219e-17,0.0024785457043284276,9.354417496250849e-17,0.0024546904085875065,1.2488589642405877e-16,0.0024428775264784264],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}},"voltage":{"waveform":{"data":[-20.5,70.5,70.5,-20.5,-20.5],"time":[0,0,0.000005,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":100,"offset":0,"label":"Rectangular","acEffectiveFrequency":591485.5360118389,"effectiveFrequency":553357.3374711702,"peak":70.5,"rms":51.572309672924284,"thd":0.4833151484524849},"harmonics":{"amplitudes":[24.2890625,57.92076613061847,1.421875,19.27588907896988,1.421875,11.528257939603122,1.421875,8.194467538528329,1.421875,6.331896912839248,1.421875,5.137996046859012,1.421875,4.304077056139349,1.421875,3.6860723299088454,1.421875,3.207698601961777,1.421875,2.8247804703632298,1.421875,2.509960393415185,1.421875,2.2453859950684323,1.421875,2.01890737840567,1.421875,1.8219644341144872,1.421875,1.6483482744897402,1.421875,1.4934420157473332,1.421875,1.3537375367153817,1.421875,1.2265178099275544,1.421875,1.1096421410704556,1.421875,1.0013973584174929,1.421875,0.9003924136274832,1.421875,0.8054822382572133,1.421875,0.7157117294021269,1.421875,0.6302738400635857,1.421875,0.5484777114167545,1.421875,0.46972405216147894,1.421875,0.3934858059809043,1.421875,0.31929270856030145,1.421875,0.24671871675852053,1.421875,0.17537155450693565,1.421875,0.10488380107099537,1.421875,0.034905072061178544],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}}},{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"ancillaryLabel":null,"data":[-5,5,-5],"numberPeriods":null,"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"Triangular","acEffectiveFrequency":110746.40291779551,"effectiveFrequency":110746.40291779551,"peak":5,"rms":2.8874560332150576,"thd":0.12151487440704967},"harmonics":{"amplitudes":[1.1608769501236793e-14,4.05366124583194,1.787369544444173e-15,0.4511310569983995,9.749053004706756e-16,0.16293015292554872,4.036157626725542e-16,0.08352979924600704,3.4998295008010614e-16,0.0508569581336163,3.1489164048780735e-16,0.034320410449418075,3.142469873118059e-16,0.024811988673843106,2.3653352035940994e-16,0.018849001010678823,2.9306524147249266e-16,0.014866633059596499,1.796485796132569e-16,0.012077180559557785,1.6247782523152451e-16,0.010049063750920609,1.5324769149805092e-16,0.008529750975091871,1.0558579733068502e-16,0.007363501410705499,7.513269775674661e-17,0.006450045785294609,5.871414177162291e-17,0.005722473794997712,9.294731722001391e-17,0.005134763398167541,1.194820309200107e-16,0.004654430423785411,8.2422739080512e-17,0.004258029771397032,9.5067306351894e-17,0.0039283108282380024,1.7540347128474968e-16,0.0036523670873925395,9.623794010508822e-17,0.0034204021424253787,1.4083470894369491e-16,0.0032248884817922927,1.4749333016985644e-16,0.0030599828465501895,1.0448590642474364e-16,0.002921112944200383,7.575487373767413e-18,0.002804680975178716,7.419510610361002e-17,0.0027078483284668376,3.924741709073613e-17,0.0026283777262804953,2.2684279102637236e-17,0.0025645167846443107,8.997077625295079e-17,0.0025149120164513483,7.131074184849219e-17,0.0024785457043284276,9.354417496250849e-17,0.0024546904085875065,1.2488589642405877e-16,0.0024428775264784264],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}},"voltage":{"waveform":{"ancillaryLabel":null,"data":[-50,50,50,-50,-50],"numberPeriods":null,"time":[0,0,0.000005,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":100,"offset":0,"label":"Rectangular","acEffectiveFrequency":591485.536011839,"effectiveFrequency":591449.4202715514,"peak":50,"rms":50,"thd":0.48331514845248497},"harmonics":{"amplitudes":[0.78125,63.64919355013018,1.5625,21.18229569117569,1.5625,12.668415318245188,1.5625,9.004909382998164,1.5625,6.958128475647527,1.5625,5.646149502042871,1.5625,4.729755006746538,1.5625,4.050628933965765,1.5625,3.524943518639316,1.5625,3.104154363036517,1.5625,2.7581982345221827,1.5625,2.467457137437843,1.5625,2.2185795367095267,1.5625,2.0021587188071255,1.5625,1.8113717302085082,1.5625,1.6411450722498175,1.5625,1.487623666720196,1.5625,1.3478217691511587,1.5625,1.2193869682092893,1.5625,1.100436657601639,1.5625,0.9894422127774558,1.5625,0.8851453167661671,1.5625,0.7864964059364037,1.5625,0.6926086154544899,1.5625,0.60272275979863,1.5625,0.5161802771005264,1.5625,0.43240198459440116,1.5625,0.3508711083080249,1.5625,0.27111946896540395,1.5625,0.192715993963664,1.5625,0.11525692425384548,1.5625,0.03835722204524927],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}}}]}]},"magnetic":{"coil":{"bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.00635,"columnShape":"rectangular","columnThickness":0,"columnWidth":0.0011999999999999997,"coordinates":[0,0,0],"pins":null,"wallThickness":0,"windingWindows":[{"angle":360,"area":0.00004901669937763475,"coordinates":[0.00395,0,0],"height":null,"radialHeight":0.00395,"sectionsAlignment":null,"sectionsOrientation":null,"shape":"round","width":null}]}},"functionalDescription":[{"isolationSide":"primary","name":"Primary","numberParallels":1,"numberTurns":45,"wire":"Dummy"},{"name":"Secondary","numberTurns":45,"numberParallels":1,"isolationSide":"secondary","wire":"Dummy"}],"layersDescription":null,"sectionsDescription":null,"turnsDescription":null},"core":{"distributorsInfo":[],"functionalDescription":{"coating":null,"gapping":[],"material":"67","numberStacks":1,"shape":{"aliases":[],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0127},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0079},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0127}},"family":"t","familySubtype":null,"magneticCircuit":"closed","name":"T 12.7/7.9/12.7","type":"standard"},"type":"toroidal","magneticCircuit":"closed"},"geometricalDescription":[{"coordinates":[0,0,0],"dimensions":null,"insulationMaterial":null,"machining":null,"material":"67","rotation":[1.5707963267948966,1.5707963267948966,0],"shape":{"aliases":[],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0127},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0079},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0127}},"family":"t","familySubtype":null,"magneticCircuit":"closed","name":"T 12.7/7.9/12.7","type":"standard"},"type":"toroidal"}],"manufacturerInfo":{"cost":null,"datasheetUrl":"https://fair-rite.com/product/toroids-5967001901/","family":null,"name":"Fair-Rite","orderCode":null,"reference":"5967001901","status":"production"},"name":"T 12.7/7.9/12.7 - 67 - Ungapped","processedDescription":{"columns":[{"area":0.000031,"coordinates":[0,0,0],"depth":0.0127,"height":0.03235840433197487,"minimumDepth":null,"minimumWidth":null,"shape":"rectangular","type":"central","width":0.0023999999999999994}],"depth":0.0127,"effectiveParameters":{"effectiveArea":0.000030479999999999996,"effectiveLength":0.03235840433197488,"effectiveVolume":9.862841640385941e-7,"minimumArea":0.00003047999999999999},"height":0.0127,"width":0.0127,"windingWindows":[{"angle":360,"area":0.00004901669937763475,"coordinates":[0.0023999999999999994,0],"height":null,"radialHeight":0.00395,"sectionsAlignment":null,"sectionsOrientation":null,"shape":null,"width":null}]}},"manufacturerInfo":{"name":"OpenMagnetics","reference":"My custom magnetic"}},"outputs":[]})"));

    for (size_t windingIndex = 0; windingIndex < mas.get_magnetic().get_coil().get_functional_description().size(); ++windingIndex) {
        mas.get_mutable_magnetic().get_mutable_coil().get_mutable_functional_description()[windingIndex].set_wire("Dummy");
    }
    mas.get_mutable_magnetic().get_mutable_coil().set_turns_description(std::nullopt);
    mas.get_mutable_magnetic().get_mutable_coil().set_layers_description(std::nullopt);
    mas.get_mutable_magnetic().get_mutable_coil().set_sections_description(std::nullopt);

    CoilAdviser coilAdviser;
    auto masMagneticsWithCoil = coilAdviser.get_advised_coil(mas, 1);
}

TEST_CASE("Test_CoilAdviser_Random_16", "[adviser][coil-adviser][bug][smoke-test]") {
    OpenMagnetics::Mas mas(json::parse(R"({"inputs":{"designRequirements":{"magnetizingInductance":{"nominal":0.0001},"minimumImpedance":[{"frequency":100000,"impedance":{"magnitude":1000}}],"name":"My Design Requirements","turnsRatios":[{"nominal":1}]},"operatingPoints":[{"name":"Operating Point No. 1","conditions":{"ambientTemperature":42},"excitationsPerWinding":[{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"data":[-5,5,-5],"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"Triangular","acEffectiveFrequency":110746.40291779551,"effectiveFrequency":110746.40291779551,"peak":5,"rms":2.8874560332150576,"thd":0.12151487440704967},"harmonics":{"amplitudes":[1.1608769501236793e-14,4.05366124583194,1.787369544444173e-15,0.4511310569983995,9.749053004706756e-16,0.16293015292554872,4.036157626725542e-16,0.08352979924600704,3.4998295008010614e-16,0.0508569581336163,3.1489164048780735e-16,0.034320410449418075,3.142469873118059e-16,0.024811988673843106,2.3653352035940994e-16,0.018849001010678823,2.9306524147249266e-16,0.014866633059596499,1.796485796132569e-16,0.012077180559557785,1.6247782523152451e-16,0.010049063750920609,1.5324769149805092e-16,0.008529750975091871,1.0558579733068502e-16,0.007363501410705499,7.513269775674661e-17,0.006450045785294609,5.871414177162291e-17,0.005722473794997712,9.294731722001391e-17,0.005134763398167541,1.194820309200107e-16,0.004654430423785411,8.2422739080512e-17,0.004258029771397032,9.5067306351894e-17,0.0039283108282380024,1.7540347128474968e-16,0.0036523670873925395,9.623794010508822e-17,0.0034204021424253787,1.4083470894369491e-16,0.0032248884817922927,1.4749333016985644e-16,0.0030599828465501895,1.0448590642474364e-16,0.002921112944200383,7.575487373767413e-18,0.002804680975178716,7.419510610361002e-17,0.0027078483284668376,3.924741709073613e-17,0.0026283777262804953,2.2684279102637236e-17,0.0025645167846443107,8.997077625295079e-17,0.0025149120164513483,7.131074184849219e-17,0.0024785457043284276,9.354417496250849e-17,0.0024546904085875065,1.2488589642405877e-16,0.0024428775264784264],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}},"voltage":{"waveform":{"data":[-20.5,70.5,70.5,-20.5,-20.5],"time":[0,0,0.000005,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":100,"offset":0,"label":"Rectangular","acEffectiveFrequency":591485.5360118389,"effectiveFrequency":553357.3374711702,"peak":70.5,"rms":51.572309672924284,"thd":0.4833151484524849},"harmonics":{"amplitudes":[24.2890625,57.92076613061847,1.421875,19.27588907896988,1.421875,11.528257939603122,1.421875,8.194467538528329,1.421875,6.331896912839248,1.421875,5.137996046859012,1.421875,4.304077056139349,1.421875,3.6860723299088454,1.421875,3.207698601961777,1.421875,2.8247804703632298,1.421875,2.509960393415185,1.421875,2.2453859950684323,1.421875,2.01890737840567,1.421875,1.8219644341144872,1.421875,1.6483482744897402,1.421875,1.4934420157473332,1.421875,1.3537375367153817,1.421875,1.2265178099275544,1.421875,1.1096421410704556,1.421875,1.0013973584174929,1.421875,0.9003924136274832,1.421875,0.8054822382572133,1.421875,0.7157117294021269,1.421875,0.6302738400635857,1.421875,0.5484777114167545,1.421875,0.46972405216147894,1.421875,0.3934858059809043,1.421875,0.31929270856030145,1.421875,0.24671871675852053,1.421875,0.17537155450693565,1.421875,0.10488380107099537,1.421875,0.034905072061178544],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}}},{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"ancillaryLabel":null,"data":[-5,5,-5],"numberPeriods":null,"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"Triangular","acEffectiveFrequency":110746.40291779551,"effectiveFrequency":110746.40291779551,"peak":5,"rms":2.8874560332150576,"thd":0.12151487440704967},"harmonics":{"amplitudes":[1.1608769501236793e-14,4.05366124583194,1.787369544444173e-15,0.4511310569983995,9.749053004706756e-16,0.16293015292554872,4.036157626725542e-16,0.08352979924600704,3.4998295008010614e-16,0.0508569581336163,3.1489164048780735e-16,0.034320410449418075,3.142469873118059e-16,0.024811988673843106,2.3653352035940994e-16,0.018849001010678823,2.9306524147249266e-16,0.014866633059596499,1.796485796132569e-16,0.012077180559557785,1.6247782523152451e-16,0.010049063750920609,1.5324769149805092e-16,0.008529750975091871,1.0558579733068502e-16,0.007363501410705499,7.513269775674661e-17,0.006450045785294609,5.871414177162291e-17,0.005722473794997712,9.294731722001391e-17,0.005134763398167541,1.194820309200107e-16,0.004654430423785411,8.2422739080512e-17,0.004258029771397032,9.5067306351894e-17,0.0039283108282380024,1.7540347128474968e-16,0.0036523670873925395,9.623794010508822e-17,0.0034204021424253787,1.4083470894369491e-16,0.0032248884817922927,1.4749333016985644e-16,0.0030599828465501895,1.0448590642474364e-16,0.002921112944200383,7.575487373767413e-18,0.002804680975178716,7.419510610361002e-17,0.0027078483284668376,3.924741709073613e-17,0.0026283777262804953,2.2684279102637236e-17,0.0025645167846443107,8.997077625295079e-17,0.0025149120164513483,7.131074184849219e-17,0.0024785457043284276,9.354417496250849e-17,0.0024546904085875065,1.2488589642405877e-16,0.0024428775264784264],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}},"voltage":{"waveform":{"ancillaryLabel":null,"data":[-50,50,50,-50,-50],"numberPeriods":null,"time":[0,0,0.000005,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":100,"offset":0,"label":"Rectangular","acEffectiveFrequency":591485.536011839,"effectiveFrequency":591449.4202715514,"peak":50,"rms":50,"thd":0.48331514845248497},"harmonics":{"amplitudes":[0.78125,63.64919355013018,1.5625,21.18229569117569,1.5625,12.668415318245188,1.5625,9.004909382998164,1.5625,6.958128475647527,1.5625,5.646149502042871,1.5625,4.729755006746538,1.5625,4.050628933965765,1.5625,3.524943518639316,1.5625,3.104154363036517,1.5625,2.7581982345221827,1.5625,2.467457137437843,1.5625,2.2185795367095267,1.5625,2.0021587188071255,1.5625,1.8113717302085082,1.5625,1.6411450722498175,1.5625,1.487623666720196,1.5625,1.3478217691511587,1.5625,1.2193869682092893,1.5625,1.100436657601639,1.5625,0.9894422127774558,1.5625,0.8851453167661671,1.5625,0.7864964059364037,1.5625,0.6926086154544899,1.5625,0.60272275979863,1.5625,0.5161802771005264,1.5625,0.43240198459440116,1.5625,0.3508711083080249,1.5625,0.27111946896540395,1.5625,0.192715993963664,1.5625,0.11525692425384548,1.5625,0.03835722204524927],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}}}]}]},"magnetic":{"coil":{"bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.00635,"columnShape":"rectangular","columnThickness":0,"columnWidth":0.0011999999999999997,"coordinates":[0,0,0],"pins":null,"wallThickness":0,"windingWindows":[{"angle":360,"area":0.00004901669937763475,"coordinates":[0.00395,0,0],"height":null,"radialHeight":0.00395,"sectionsAlignment":null,"sectionsOrientation":null,"shape":"round","width":null}]}},"functionalDescription":[{"isolationSide":"primary","name":"Primary","numberParallels":1,"numberTurns":45,"wire":"Dummy"},{"name":"Secondary","numberTurns":45,"numberParallels":1,"isolationSide":"secondary","wire":"Dummy"}],"layersDescription":null,"sectionsDescription":null,"turnsDescription":null},"core":{"distributorsInfo":[],"functionalDescription":{"coating":null,"gapping":[],"material":"67","numberStacks":1,"shape":{"aliases":[],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0127},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0079},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0127}},"family":"t","familySubtype":null,"magneticCircuit":"closed","name":"T 12.7/7.9/12.7","type":"standard"},"type":"toroidal","magneticCircuit":"closed"},"geometricalDescription":[{"coordinates":[0,0,0],"dimensions":null,"insulationMaterial":null,"machining":null,"material":"67","rotation":[1.5707963267948966,1.5707963267948966,0],"shape":{"aliases":[],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0127},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0079},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0127}},"family":"t","familySubtype":null,"magneticCircuit":"closed","name":"T 12.7/7.9/12.7","type":"standard"},"type":"toroidal"}],"manufacturerInfo":{"cost":null,"datasheetUrl":"https://fair-rite.com/product/toroids-5967001901/","family":null,"name":"Fair-Rite","orderCode":null,"reference":"5967001901","status":"production"},"name":"T 12.7/7.9/12.7 - 67 - Ungapped","processedDescription":{"columns":[{"area":0.000031,"coordinates":[0,0,0],"depth":0.0127,"height":0.03235840433197487,"minimumDepth":null,"minimumWidth":null,"shape":"rectangular","type":"central","width":0.0023999999999999994}],"depth":0.0127,"effectiveParameters":{"effectiveArea":0.000030479999999999996,"effectiveLength":0.03235840433197488,"effectiveVolume":9.862841640385941e-7,"minimumArea":0.00003047999999999999},"height":0.0127,"width":0.0127,"windingWindows":[{"angle":360,"area":0.00004901669937763475,"coordinates":[0.0023999999999999994,0],"height":null,"radialHeight":0.00395,"sectionsAlignment":null,"sectionsOrientation":null,"shape":null,"width":null}]}},"manufacturerInfo":{"name":"OpenMagnetics","reference":"My custom magnetic"}},"outputs":[]})"));

    for (size_t windingIndex = 0; windingIndex < mas.get_magnetic().get_coil().get_functional_description().size(); ++windingIndex) {
        mas.get_mutable_magnetic().get_mutable_coil().get_mutable_functional_description()[windingIndex].set_wire("Dummy");
    }
    mas.get_mutable_magnetic().get_mutable_coil().set_turns_description(std::nullopt);
    mas.get_mutable_magnetic().get_mutable_coil().set_layers_description(std::nullopt);
    mas.get_mutable_magnetic().get_mutable_coil().set_sections_description(std::nullopt);

    CoilAdviser coilAdviser;
    auto masMagneticsWithCoil = coilAdviser.get_advised_coil(mas, 1);
}

TEST_CASE("Test_CoilAdviser_Random_17", "[adviser][coil-adviser][bug][smoke-test]") {
    auto json_path_1964 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_coiladviser_random_17_1964.json");
    std::ifstream json_file_1964(json_path_1964);
    OpenMagnetics::Mas mas(json::parse(json_file_1964));

    for (size_t windingIndex = 0; windingIndex < mas.get_magnetic().get_coil().get_functional_description().size(); ++windingIndex) {
        mas.get_mutable_magnetic().get_mutable_coil().get_mutable_functional_description()[windingIndex].set_wire("Dummy");
    }
    mas.get_mutable_magnetic().get_mutable_coil().set_turns_description(std::nullopt);
    mas.get_mutable_magnetic().get_mutable_coil().set_layers_description(std::nullopt);
    mas.get_mutable_magnetic().get_mutable_coil().set_sections_description(std::nullopt);

    CoilAdviser coilAdviser;
    auto masMagneticsWithCoil = coilAdviser.get_advised_coil(mas, 1);

    json result = json();
    to_json(result, masMagneticsWithCoil[0]);

    auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
    auto outFile = outputFilePath;
    std::string filename = "Test_CoilAdviser_Random_17.svg";
    outFile.append(filename);
    Painter painter(outFile);

    painter.paint_core(masMagneticsWithCoil[0].get_mutable_magnetic());
    painter.paint_bobbin(masMagneticsWithCoil[0].get_mutable_magnetic());
    painter.paint_coil_turns(masMagneticsWithCoil[0].get_mutable_magnetic());
    painter.export_svg();
}

TEST_CASE("Test_SectionsAdviser_Web_0", "[adviser][coil-adviser][bug][smoke-test]") {
    std::string patternString = "[0]";
    json patternJson = json::parse(patternString);
    std::vector<size_t> pattern; 
    for (auto& elem : patternJson) {
        pattern.push_back(elem);
    }
    int repetitions = 1;
    OpenMagnetics::Mas mas(json::parse(R"({"inputs":{"designRequirements":{"magnetizingInductance":{"nominal":0.0001},"name":"My Design Requirements","turnsRatios":[]},"operatingPoints":[{"name":"Operating Point No. 1","conditions":{"ambientTemperature":42},"excitationsPerWinding":[{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"ancillaryLabel":null,"data":[-5,5,-5],"numberPeriods":null,"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"Triangular","acEffectiveFrequency":110746.40291779551,"effectiveFrequency":110746.40291779551,"peak":5,"rms":2.8874560332150576,"thd":0.12151487440704967},"harmonics":{"amplitudes":[1.1608769501236793e-14,4.05366124583194,1.787369544444173e-15,0.4511310569983995,9.749053004706756e-16,0.16293015292554872,4.036157626725542e-16,0.08352979924600704,3.4998295008010614e-16,0.0508569581336163,3.1489164048780735e-16,0.034320410449418075,3.142469873118059e-16,0.024811988673843106,2.3653352035940994e-16,0.018849001010678823,2.9306524147249266e-16,0.014866633059596499,1.796485796132569e-16,0.012077180559557785,1.6247782523152451e-16,0.010049063750920609,1.5324769149805092e-16,0.008529750975091871,1.0558579733068502e-16,0.007363501410705499,7.513269775674661e-17,0.006450045785294609,5.871414177162291e-17,0.005722473794997712,9.294731722001391e-17,0.005134763398167541,1.194820309200107e-16,0.004654430423785411,8.2422739080512e-17,0.004258029771397032,9.5067306351894e-17,0.0039283108282380024,1.7540347128474968e-16,0.0036523670873925395,9.623794010508822e-17,0.0034204021424253787,1.4083470894369491e-16,0.0032248884817922927,1.4749333016985644e-16,0.0030599828465501895,1.0448590642474364e-16,0.002921112944200383,7.575487373767413e-18,0.002804680975178716,7.419510610361002e-17,0.0027078483284668376,3.924741709073613e-17,0.0026283777262804953,2.2684279102637236e-17,0.0025645167846443107,8.997077625295079e-17,0.0025149120164513483,7.131074184849219e-17,0.0024785457043284276,9.354417496250849e-17,0.0024546904085875065,1.2488589642405877e-16,0.0024428775264784264],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}},"voltage":{"waveform":{"ancillaryLabel":null,"data":[-50,50,50,-50,-50],"numberPeriods":null,"time":[0,0,0.000005,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":100,"offset":0,"label":"Rectangular","acEffectiveFrequency":591485.536011839,"effectiveFrequency":591449.4202715514,"peak":50,"rms":50,"thd":0.48331514845248497},"harmonics":{"amplitudes":[0.78125,63.64919355013018,1.5625,21.18229569117569,1.5625,12.668415318245188,1.5625,9.004909382998164,1.5625,6.958128475647527,1.5625,5.646149502042871,1.5625,4.729755006746538,1.5625,4.050628933965765,1.5625,3.524943518639316,1.5625,3.104154363036517,1.5625,2.7581982345221827,1.5625,2.467457137437843,1.5625,2.2185795367095267,1.5625,2.0021587188071255,1.5625,1.8113717302085082,1.5625,1.6411450722498175,1.5625,1.487623666720196,1.5625,1.3478217691511587,1.5625,1.2193869682092893,1.5625,1.100436657601639,1.5625,0.9894422127774558,1.5625,0.8851453167661671,1.5625,0.7864964059364037,1.5625,0.6926086154544899,1.5625,0.60272275979863,1.5625,0.5161802771005264,1.5625,0.43240198459440116,1.5625,0.3508711083080249,1.5625,0.27111946896540395,1.5625,0.192715993963664,1.5625,0.11525692425384548,1.5625,0.03835722204524927],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}}},{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"data":[-5,5,-5],"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"Triangular"}},"voltage":{"waveform":{"data":[-20.5,70.5,70.5,-20.5,-20.5],"time":[0,0,0.000005,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":100,"offset":0,"label":"Rectangular"}}},{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"data":[-5,5,-5],"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"Triangular"}},"voltage":{"waveform":{"data":[-20.5,70.5,70.5,-20.5,-20.5],"time":[0,0,0.000005,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":100,"offset":0,"label":"Rectangular"}}},{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"data":[-5,5,-5],"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"Triangular"}},"voltage":{"waveform":{"data":[-20.5,70.5,70.5,-20.5,-20.5],"time":[0,0,0.000005,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":100,"offset":0,"label":"Rectangular"}}},{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"data":[-5,5,-5],"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"Triangular"}},"voltage":{"waveform":{"data":[-20.5,70.5,70.5,-20.5,-20.5],"time":[0,0,0.000005,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":100,"offset":0,"label":"Rectangular"}}},{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"data":[-5,5,-5],"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"Triangular"}},"voltage":{"waveform":{"data":[-20.5,70.5,70.5,-20.5,-20.5],"time":[0,0,0.000005,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":100,"offset":0,"label":"Rectangular"}}},{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"data":[-5,5,-5],"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"Triangular"}},"voltage":{"waveform":{"data":[-20.5,70.5,70.5,-20.5,-20.5],"time":[0,0,0.000005,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":100,"offset":0,"label":"Rectangular"}}}]}]},"magnetic":{"coil":{"bobbin":"Dummy","functionalDescription":[{"connections":null,"isolationSide":"primary","name":"Primary","numberParallels":1,"numberTurns":1,"wire":"Dummy"}]},"core":{"distributorsInfo":null,"functionalDescription":{"coating":null,"gapping":[{"length":0.001,"type":"subtractive"},{"length":0.001,"type":"subtractive"},{"length":0.001,"type":"subtractive"},{"length":0.000005,"type":"residual"},{"length":0.000005,"type":"residual"}],"material":"3C98","numberStacks":1,"shape":{"aliases":[],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0414,"minimum":0.0396,"nominal":0.0405},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.01525,"minimum":0.01505,"nominal":0.01515},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0286,"minimum":0.0274,"nominal":0.028},"D":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.01015,"minimum":0.00985,"nominal":0.01},"E":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0376,"minimum":0.0364,"nominal":0.037},"F":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0152,"minimum":0.0146,"nominal":0.0149},"G":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":0.028,"nominal":null}},"family":"pq","familySubtype":null,"magneticCircuit":"open","name":"PQ 40/30","type":"standard"},"type":"two-piece set"},"geometricalDescription":[{"coordinates":[0,0,0],"dimensions":null,"insulationMaterial":null,"machining":null,"material":"","rotation":[3.141592653589793,3.141592653589793,0],"shape":{"aliases":[],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0414,"minimum":0.0396,"nominal":0.0405},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.01525,"minimum":0.01505,"nominal":0.01515},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0286,"minimum":0.0274,"nominal":0.028},"D":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.01015,"minimum":0.00985,"nominal":0.01},"E":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0376,"minimum":0.0364,"nominal":0.037},"F":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0152,"minimum":0.0146,"nominal":0.0149},"G":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":0.028,"nominal":null}},"family":"pq","familySubtype":null,"magneticCircuit":"open","name":"PQ 40/30","type":"standard"},"type":"half set"},{"coordinates":[0,0,0],"dimensions":null,"insulationMaterial":null,"machining":null,"material":"","rotation":[0,0,0],"shape":{"aliases":[],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0414,"minimum":0.0396,"nominal":0.0405},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.01525,"minimum":0.01505,"nominal":0.01515},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0286,"minimum":0.0274,"nominal":0.028},"D":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.01015,"minimum":0.00985,"nominal":0.01},"E":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0376,"minimum":0.0364,"nominal":0.037},"F":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0152,"minimum":0.0146,"nominal":0.0149},"G":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":0.028,"nominal":null}},"family":"pq","familySubtype":null,"magneticCircuit":"open","name":"PQ 40/30","type":"standard"},"type":"half set"}],"manufacturerInfo":null,"name":"DummyCore","processedDescription":{"columns":[{"area":0.000175,"coordinates":[0,0,0],"depth":0.0149,"height":0.02,"minimumDepth":null,"minimumWidth":null,"shape":"round","type":"central","width":0.0149},{"area":0.000101,"coordinates":[0.020304,0,0],"depth":0.028,"height":0.02,"minimumDepth":null,"minimumWidth":0.001751,"shape":"irregular","type":"lateral","width":0.003608},{"area":0.000101,"coordinates":[-0.020304,0,0],"depth":0.028,"height":0.02,"minimumDepth":null,"minimumWidth":0.001751,"shape":"irregular","type":"lateral","width":0.003608}],"depth":0.028,"effectiveParameters":{"effectiveArea":0.00019510490157384353,"effectiveLength":0.0740608824164849,"effectiveVolume":0.000014449641174340287,"minimumArea":0.0001743662462558675},"height":0.0303,"width":0.0405,"windingWindows":[{"angle":null,"area":0.00022099999999999998,"coordinates":[0.00745,0],"height":0.02,"radialHeight":null,"sectionsOrientation":null,"shape":null,"width":0.011049999999999999}]}}},"outputs":[]})"));
    mas.get_mutable_magnetic().get_mutable_coil().set_bobbin(OpenMagnetics::Bobbin::create_quick_bobbin(mas.get_mutable_magnetic().get_mutable_core()));
    auto sections = CoilAdviser().get_advised_sections(mas, pattern, repetitions);
}

TEST_CASE("Test_SectionsAdviser_Web_1", "[adviser][coil-adviser][bug][smoke-test]") {
    std::string patternString = "[0, 1]";
    json patternJson = json::parse(patternString);
    std::vector<size_t> pattern; 
    for (auto& elem : patternJson) {
        pattern.push_back(elem);
    }
    int repetitions = 1;
    OpenMagnetics::Mas mas(json::parse(R"({"inputs":{"designRequirements":{"insulation":{"altitude":{"maximum":2000},"cti":"Group IIIB","pollutionDegree":"P2","overvoltageCategory":"OVC-III","insulationType":"Double","mainSupplyVoltage":{"maximum":400},"standards":["IEC 60664-1"]},"magnetizingInductance":{"nominal":0.0001},"name":"My Design Requirements","turnsRatios":[{"nominal":1}]},"operatingPoints":[{"name":"Operating Point No. 1","conditions":{"ambientTemperature":42},"excitationsPerWinding":[{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"data":[-5,5,-5],"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"Triangular","acEffectiveFrequency":110746.40291779551,"effectiveFrequency":110746.40291779551,"peak":5,"rms":2.8874560332150576,"thd":0.12151487440704967},"harmonics":{"amplitudes":[1.1608769501236793e-14,4.05366124583194,1.787369544444173e-15,0.4511310569983995,9.749053004706756e-16,0.16293015292554872,4.036157626725542e-16,0.08352979924600704,3.4998295008010614e-16,0.0508569581336163,3.1489164048780735e-16,0.034320410449418075,3.142469873118059e-16,0.024811988673843106,2.3653352035940994e-16,0.018849001010678823,2.9306524147249266e-16,0.014866633059596499,1.796485796132569e-16,0.012077180559557785,1.6247782523152451e-16,0.010049063750920609,1.5324769149805092e-16,0.008529750975091871,1.0558579733068502e-16,0.007363501410705499,7.513269775674661e-17,0.006450045785294609,5.871414177162291e-17,0.005722473794997712,9.294731722001391e-17,0.005134763398167541,1.194820309200107e-16,0.004654430423785411,8.2422739080512e-17,0.004258029771397032,9.5067306351894e-17,0.0039283108282380024,1.7540347128474968e-16,0.0036523670873925395,9.623794010508822e-17,0.0034204021424253787,1.4083470894369491e-16,0.0032248884817922927,1.4749333016985644e-16,0.0030599828465501895,1.0448590642474364e-16,0.002921112944200383,7.575487373767413e-18,0.002804680975178716,7.419510610361002e-17,0.0027078483284668376,3.924741709073613e-17,0.0026283777262804953,2.2684279102637236e-17,0.0025645167846443107,8.997077625295079e-17,0.0025149120164513483,7.131074184849219e-17,0.0024785457043284276,9.354417496250849e-17,0.0024546904085875065,1.2488589642405877e-16,0.0024428775264784264],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}},"voltage":{"waveform":{"data":[-20.5,70.5,70.5,-20.5,-20.5],"time":[0,0,0.000005,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":100,"offset":0,"label":"Rectangular","acEffectiveFrequency":591485.5360118389,"effectiveFrequency":553357.3374711702,"peak":70.5,"rms":51.572309672924284,"thd":0.4833151484524849},"harmonics":{"amplitudes":[24.2890625,57.92076613061847,1.421875,19.27588907896988,1.421875,11.528257939603122,1.421875,8.194467538528329,1.421875,6.331896912839248,1.421875,5.137996046859012,1.421875,4.304077056139349,1.421875,3.6860723299088454,1.421875,3.207698601961777,1.421875,2.8247804703632298,1.421875,2.509960393415185,1.421875,2.2453859950684323,1.421875,2.01890737840567,1.421875,1.8219644341144872,1.421875,1.6483482744897402,1.421875,1.4934420157473332,1.421875,1.3537375367153817,1.421875,1.2265178099275544,1.421875,1.1096421410704556,1.421875,1.0013973584174929,1.421875,0.9003924136274832,1.421875,0.8054822382572133,1.421875,0.7157117294021269,1.421875,0.6302738400635857,1.421875,0.5484777114167545,1.421875,0.46972405216147894,1.421875,0.3934858059809043,1.421875,0.31929270856030145,1.421875,0.24671871675852053,1.421875,0.17537155450693565,1.421875,0.10488380107099537,1.421875,0.034905072061178544],"frequencies":[0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000,1700000,1800000,1900000,2000000,2100000,2200000,2300000,2400000,2500000,2600000,2700000,2800000,2900000,3000000,3100000,3200000,3300000,3400000,3500000,3600000,3700000,3800000,3900000,4000000,4100000,4200000,4300000,4400000,4500000,4600000,4700000,4800000,4900000,5000000,5100000,5200000,5300000,5400000,5500000,5600000,5700000,5800000,5900000,6000000,6100000,6200000,6300000]}}},{"name":"Primary winding excitation","frequency":100000,"current":{"waveform":{"data":[-5,5,-5],"time":[0,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":10,"offset":0,"label":"Triangular"}},"voltage":{"waveform":{"data":[-20.5,70.5,70.5,-20.5,-20.5],"time":[0,0,0.000005,0.000005,0.00001]},"processed":{"dutyCycle":0.5,"peakToPeak":100,"offset":0,"label":"Rectangular"}}}]}]},"magnetic":{"coil":{"bobbin":"Dummy","functionalDescription":[{"isolationSide":"primary","name":"Primary","numberParallels":1,"numberTurns":24,"wire":{"coating":{"breakdownVoltage":70,"grade":1,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":null,"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.00001},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Elektrisola","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 0.01 - Grade 1","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.000013000000000000001,"minimum":0.000012,"nominal":null},"outerHeight":null,"outerWidth":null,"standard":"IEC 60317","standardName":"0.01 mm","strand":null,"type":"round"}},{"name":"Secondary","numberTurns":1,"numberParallels":1,"isolationSide":"secondary","wire":"Dummy"}]},"core":{"distributorsInfo":[],"functionalDescription":{"coating":null,"gapping":[{"area":0.000071,"coordinates":[0,0.0005,0],"distanceClosestNormalSurface":0.01125,"distanceClosestParallelSurface":0.006625,"length":0.001,"sectionDimensions":[0.0095,0.0095],"shape":"round","type":"subtractive"},{"area":0.00006,"coordinates":[0.014533,0,0],"distanceClosestNormalSurface":0.012248,"distanceClosestParallelSurface":0.006625,"length":0.000005,"sectionDimensions":[0.006316,0.0095],"shape":"irregular","type":"residual"},{"area":0.00006,"coordinates":[-0.014533,0,0],"distanceClosestNormalSurface":0.012248,"distanceClosestParallelSurface":0.006625,"length":0.000005,"sectionDimensions":[0.006316,0.0095],"shape":"irregular","type":"residual"}],"material":"3C91","numberStacks":1,"shape":{"aliases":["EC 35/17/10"],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0353,"minimum":0.0337,"nominal":null},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.01745,"minimum":0.01715,"nominal":null},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0098,"minimum":0.0092,"nominal":null},"D":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0126,"minimum":0.0119,"nominal":null},"E":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0233,"minimum":0.0222,"nominal":null},"F":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0098,"minimum":0.0092,"nominal":null},"T":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0293,"minimum":0.0277,"nominal":null},"r":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0005},"s":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.003,"minimum":0.0025,"nominal":null}},"family":"ec","familySubtype":null,"magneticCircuit":"open","name":"EC 35/17/10","type":"standard"},"type":"two-piece set"},"geometricalDescription":[{"coordinates":[0,0,0],"dimensions":null,"insulationMaterial":null,"machining":[{"coordinates":[0,0.0005,0],"length":0.001}],"material":"3C91","rotation":[3.141592653589793,3.141592653589793,0],"shape":{"aliases":["EC 35/17/10"],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0353,"minimum":0.0337,"nominal":null},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.01745,"minimum":0.01715,"nominal":null},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0098,"minimum":0.0092,"nominal":null},"D":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0126,"minimum":0.0119,"nominal":null},"E":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0233,"minimum":0.0222,"nominal":null},"F":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0098,"minimum":0.0092,"nominal":null},"T":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0293,"minimum":0.0277,"nominal":null},"r":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0005},"s":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.003,"minimum":0.0025,"nominal":null}},"family":"ec","familySubtype":null,"magneticCircuit":"open","name":"EC 35/17/10","type":"standard"},"type":"half set"},{"coordinates":[0,0,0],"dimensions":null,"insulationMaterial":null,"machining":null,"material":"3C91","rotation":[0,0,0],"shape":{"aliases":["EC 35/17/10"],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0353,"minimum":0.0337,"nominal":null},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.01745,"minimum":0.01715,"nominal":null},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0098,"minimum":0.0092,"nominal":null},"D":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0126,"minimum":0.0119,"nominal":null},"E":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0233,"minimum":0.0222,"nominal":null},"F":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0098,"minimum":0.0092,"nominal":null},"T":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0293,"minimum":0.0277,"nominal":null},"r":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0005},"s":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.003,"minimum":0.0025,"nominal":null}},"family":"ec","familySubtype":null,"magneticCircuit":"open","name":"EC 35/17/10","type":"standard"},"type":"half set"}],"manufacturerInfo":{"cost":null,"datasheetUrl":"https://ferroxcube.com/upload/media/product/file/Pr_ds/EC35_17_10.pdf","family":null,"name":"Ferroxcube","orderCode":null,"reference":"EC35/17/10-3C91-G1000","status":"production"},"name":"EC 35/17/10 - 3C91 - Gapped 1.000 mm","processedDescription":{"columns":[{"area":0.000071,"coordinates":[0,0,0],"depth":0.0095,"height":0.0245,"minimumDepth":null,"minimumWidth":null,"shape":"round","type":"central","width":0.0095},{"area":0.00006,"coordinates":[0.014533,0,0],"depth":0.0095,"height":0.0245,"minimumDepth":null,"minimumWidth":0.005876,"shape":"irregular","type":"lateral","width":0.006316},{"area":0.00006,"coordinates":[-0.014533,0,0],"depth":0.0095,"height":0.0245,"minimumDepth":null,"minimumWidth":0.005876,"shape":"irregular","type":"lateral","width":0.006316}],"depth":0.0095,"effectiveParameters":{"effectiveArea":0.00008700294107630485,"effectiveLength":0.07610520812632617,"effectiveVolume":0.000006621376938214672,"minimumArea":0.0000708821842466197},"height":0.0346,"width":0.0345,"windingWindows":[{"angle":null,"area":0.0001623125,"coordinates":[0.00475,0],"height":0.0245,"radialHeight":null,"sectionsOrientation":null,"shape":null,"width":0.006625}]}}},"outputs":[]})"));

    auto bobbin = mas.get_magnetic().get_coil().get_bobbin();
    if (std::holds_alternative<std::string>(bobbin)) {
        auto bobbinString = std::get<std::string>(bobbin);
        if (bobbinString == "Dummy") {
            mas.get_mutable_magnetic().get_mutable_coil().set_bobbin(OpenMagnetics::Bobbin::create_quick_bobbin(mas.get_mutable_magnetic().get_mutable_core()));
        }
    }
    mas.get_mutable_inputs().process();

    auto sections = CoilAdviser().get_advised_sections(mas, pattern, repetitions);
}

TEST_CASE("Test_SectionsAdviser_Web_2", "[adviser][coil-adviser][bug][smoke-test]") {
    std::string patternString = "[0]";
    json patternJson = json::parse(patternString);
    std::vector<size_t> pattern; 
    for (auto& elem : patternJson) {
        pattern.push_back(elem);
    }
    int repetitions = 1;
    OpenMagnetics::Mas mas(json::parse(R"({"inputs":{"designRequirements":{"insulation":null,"isolationSides":null,"leakageInductance":null,"magnetizingInductance":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":0.0001},"market":null,"maximumDimensions":null,"maximumWeight":null,"name":"My Design Requirements","operatingTemperature":null,"strayCapacitance":null,"terminalType":null,"topology":null,"turnsRatios":[],"wiringTechnology":null},"operatingPoints":[{"conditions":{"ambientRelativeHumidity":null,"ambientTemperature":42.0,"cooling":null,"name":null},"excitationsPerWinding":[{"current":{"harmonics":{"amplitudes":[1.1608769501236793e-14,4.05366124583194,1.787369544444173e-15,0.4511310569983995,9.749053004706756e-16,0.16293015292554872,4.036157626725542e-16,0.08352979924600704,3.4998295008010614e-16,0.0508569581336163,3.1489164048780735e-16,0.034320410449418075,3.142469873118059e-16,0.024811988673843106,2.3653352035940994e-16,0.018849001010678823,2.9306524147249266e-16,0.014866633059596499,1.796485796132569e-16,0.012077180559557785,1.6247782523152451e-16,0.010049063750920609,1.5324769149805092e-16,0.008529750975091871,1.0558579733068502e-16,0.007363501410705499,7.513269775674661e-17,0.006450045785294609,5.871414177162291e-17,0.005722473794997712,9.294731722001391e-17,0.005134763398167541,1.194820309200107e-16,0.004654430423785411,8.2422739080512e-17,0.004258029771397032,9.5067306351894e-17,0.0039283108282380024,1.7540347128474968e-16,0.0036523670873925395,9.623794010508822e-17,0.0034204021424253787,1.4083470894369491e-16,0.0032248884817922927,1.4749333016985644e-16,0.0030599828465501895,1.0448590642474364e-16,0.002921112944200383,7.575487373767413e-18,0.002804680975178716,7.419510610361002e-17,0.0027078483284668376,3.924741709073613e-17,0.0026283777262804953,2.2684279102637236e-17,0.0025645167846443107,8.997077625295079e-17,0.0025149120164513483,7.131074184849219e-17,0.0024785457043284276,9.354417496250849e-17,0.0024546904085875065,1.2488589642405877e-16,0.0024428775264784264],"frequencies":[0.0,100000.0,200000.0,300000.0,400000.0,500000.0,600000.0,700000.0,800000.0,900000.0,1000000.0,1100000.0,1200000.0,1300000.0,1400000.0,1500000.0,1600000.0,1700000.0,1800000.0,1900000.0,2000000.0,2100000.0,2200000.0,2300000.0,2400000.0,2500000.0,2600000.0,2700000.0,2800000.0,2900000.0,3000000.0,3100000.0,3200000.0,3300000.0,3400000.0,3500000.0,3600000.0,3700000.0,3800000.0,3900000.0,4000000.0,4100000.0,4200000.0,4300000.0,4400000.0,4500000.0,4600000.0,4700000.0,4800000.0,4900000.0,5000000.0,5100000.0,5200000.0,5300000.0,5400000.0,5500000.0,5600000.0,5700000.0,5800000.0,5900000.0,6000000.0,6100000.0,6200000.0,6300000.0]},"processed":{"acEffectiveFrequency":110746.40291779551,"average":null,"dutyCycle":0.5,"effectiveFrequency":110746.40291779551,"label":"Triangular","offset":0.0,"peak":5.0,"peakToPeak":10.0,"phase":null,"rms":2.8874560332150576,"thd":0.12151487440704967},"waveform":{"ancillaryLabel":null,"data":[-5.0,5.0,-5.0],"numberPeriods":null,"time":[0.0,5e-06,1e-05]}},"frequency":100000.0,"magneticFieldStrength":null,"magneticFluxDensity":null,"magnetizingCurrent":null,"name":"Primary winding excitation","voltage":{"harmonics":{"amplitudes":[24.2890625,57.92076613061847,1.421875,19.27588907896988,1.421875,11.528257939603122,1.421875,8.194467538528329,1.421875,6.331896912839248,1.421875,5.137996046859012,1.421875,4.304077056139349,1.421875,3.6860723299088454,1.421875,3.207698601961777,1.421875,2.8247804703632298,1.421875,2.509960393415185,1.421875,2.2453859950684323,1.421875,2.01890737840567,1.421875,1.8219644341144872,1.421875,1.6483482744897402,1.421875,1.4934420157473332,1.421875,1.3537375367153817,1.421875,1.2265178099275544,1.421875,1.1096421410704556,1.421875,1.0013973584174929,1.421875,0.9003924136274832,1.421875,0.8054822382572133,1.421875,0.7157117294021269,1.421875,0.6302738400635857,1.421875,0.5484777114167545,1.421875,0.46972405216147894,1.421875,0.3934858059809043,1.421875,0.31929270856030145,1.421875,0.24671871675852053,1.421875,0.17537155450693565,1.421875,0.10488380107099537,1.421875,0.034905072061178544],"frequencies":[0.0,100000.0,200000.0,300000.0,400000.0,500000.0,600000.0,700000.0,800000.0,900000.0,1000000.0,1100000.0,1200000.0,1300000.0,1400000.0,1500000.0,1600000.0,1700000.0,1800000.0,1900000.0,2000000.0,2100000.0,2200000.0,2300000.0,2400000.0,2500000.0,2600000.0,2700000.0,2800000.0,2900000.0,3000000.0,3100000.0,3200000.0,3300000.0,3400000.0,3500000.0,3600000.0,3700000.0,3800000.0,3900000.0,4000000.0,4100000.0,4200000.0,4300000.0,4400000.0,4500000.0,4600000.0,4700000.0,4800000.0,4900000.0,5000000.0,5100000.0,5200000.0,5300000.0,5400000.0,5500000.0,5600000.0,5700000.0,5800000.0,5900000.0,6000000.0,6100000.0,6200000.0,6300000.0]},"processed":{"acEffectiveFrequency":591485.5360118389,"average":null,"dutyCycle":0.5,"effectiveFrequency":553357.3374711702,"label":"Rectangular","offset":0.0,"peak":70.5,"peakToPeak":100.0,"phase":null,"rms":51.572309672924284,"thd":0.4833151484524849},"waveform":{"ancillaryLabel":null,"data":[-20.5,70.5,70.5,-20.5,-20.5],"numberPeriods":null,"time":[0.0,0.0,5e-06,5e-06,1e-05]}}},{"current":{"harmonics":null,"processed":{"acEffectiveFrequency":null,"average":null,"dutyCycle":0.5,"effectiveFrequency":null,"label":"Triangular","offset":0.0,"peak":null,"peakToPeak":10.0,"phase":null,"rms":null,"thd":null},"waveform":{"ancillaryLabel":null,"data":[-5.0,5.0,-5.0],"numberPeriods":null,"time":[0.0,5e-06,1e-05]}},"frequency":100000.0,"magneticFieldStrength":null,"magneticFluxDensity":null,"magnetizingCurrent":null,"name":"Primary winding excitation","voltage":{"harmonics":null,"processed":{"acEffectiveFrequency":null,"average":null,"dutyCycle":0.5,"effectiveFrequency":null,"label":"Rectangular","offset":0.0,"peak":null,"peakToPeak":100.0,"phase":null,"rms":null,"thd":null},"waveform":{"ancillaryLabel":null,"data":[-20.5,70.5,70.5,-20.5,-20.5],"numberPeriods":null,"time":[0.0,0.0,5e-06,5e-06,1e-05]}}},{"current":{"harmonics":null,"processed":{"acEffectiveFrequency":null,"average":null,"dutyCycle":0.5,"effectiveFrequency":null,"label":"Triangular","offset":0.0,"peak":null,"peakToPeak":10.0,"phase":null,"rms":null,"thd":null},"waveform":{"ancillaryLabel":null,"data":[-5.0,5.0,-5.0],"numberPeriods":null,"time":[0.0,5e-06,1e-05]}},"frequency":100000.0,"magneticFieldStrength":null,"magneticFluxDensity":null,"magnetizingCurrent":null,"name":"Primary winding excitation","voltage":{"harmonics":null,"processed":{"acEffectiveFrequency":null,"average":null,"dutyCycle":0.5,"effectiveFrequency":null,"label":"Rectangular","offset":0.0,"peak":null,"peakToPeak":100.0,"phase":null,"rms":null,"thd":null},"waveform":{"ancillaryLabel":null,"data":[-20.5,70.5,70.5,-20.5,-20.5],"numberPeriods":null,"time":[0.0,0.0,5e-06,5e-06,1e-05]}}}],"name":"Operating Point No. 1"}]},"magnetic":{"coil":{"bobbin":{"distributorsInfo":null,"functionalDescription":null,"manufacturerInfo":null,"name":null,"processedDescription":{"columnDepth":0.003912540921738127,"columnShape":"round","columnThickness":0.0010750409217381266,"columnWidth":0.003912540921738127,"coordinates":[0.0,0.0,0.0],"pins":null,"wallThickness":0.0007636652757989004,"windingWindows":[{"angle":null,"area":2.0400047558919626e-05,"coordinates":[0.004956270460869064,0.0,0.0],"height":0.009772669448402199,"radialHeight":null,"sectionsOrientation":null,"shape":"rectangular","width":0.002087459078261874}]}},"functionalDescription":[{"connections":null,"isolationSide":"primary","name":"Primary","numberParallels":1,"numberTurns":32,"wire":{"coating":{"breakdownVoltage":70.0,"grade":1,"material":null,"numberLayers":null,"temperatureRating":null,"thickness":null,"thicknessLayers":null,"type":"enamelled"},"conductingArea":null,"conductingDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":1e-05},"conductingHeight":null,"conductingWidth":null,"edgeRadius":null,"manufacturerInfo":{"cost":null,"datasheetUrl":null,"family":null,"name":"Elektrisola","orderCode":null,"reference":null,"status":null},"material":"copper","name":"Round 0.01 - Grade 1","numberConductors":1,"outerDiameter":{"excludeMaximum":null,"excludeMinimum":null,"maximum":1.3000000000000001e-05,"minimum":1.2e-05,"nominal":null},"outerHeight":null,"outerWidth":null,"standard":"IEC 60317","standardName":"0.01 mm","strand":null,"type":"round"}}],"layersDescription":null,"sectionsDescription":null,"turnsDescription":null},"core":{"distributorsInfo":[],"functionalDescription":{"coating":null,"gapping":[{"area":2.6e-05,"coordinates":[0.0,0.0,0.0],"distanceClosestNormalSurface":0.005443,"distanceClosestParallelSurface":0.0031625000000000004,"length":0.000414,"sectionDimensions":[0.005675,0.005675],"shape":"round","type":"subtractive"},{"area":0.000163,"coordinates":[0.0,0.0,-0.007038],"distanceClosestNormalSurface":0.005648,"distanceClosestParallelSurface":0.0031625000000000004,"length":5e-06,"sectionDimensions":[0.078555,0.002075],"shape":"irregular","type":"residual"}],"material":"3C91","numberStacks":1,"shape":{"aliases":[],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0184,"minimum":0.0176,"nominal":null},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0085,"minimum":0.0083,"nominal":null},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.01125,"minimum":0.01075,"nominal":null},"D":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0058,"minimum":0.0055,"nominal":null},"E":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0124,"minimum":0.0116,"nominal":null},"F":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.00585,"minimum":0.0055,"nominal":null},"K":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.002925,"minimum":null,"nominal":null}},"family":"ep","familySubtype":null,"magneticCircuit":"open","name":"EP 17","type":"standard"},"type":"two-piece set"},"geometricalDescription":[{"coordinates":[0.0,0.0,0.0],"dimensions":null,"insulationMaterial":null,"machining":[{"coordinates":[0.0,0.0001035,0.0],"length":0.000207}],"material":"3C91","rotation":[3.141592653589793,3.141592653589793,0.0],"shape":{"aliases":[],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0184,"minimum":0.0176,"nominal":null},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0085,"minimum":0.0083,"nominal":null},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.01125,"minimum":0.01075,"nominal":null},"D":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0058,"minimum":0.0055,"nominal":null},"E":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0124,"minimum":0.0116,"nominal":null},"F":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.00585,"minimum":0.0055,"nominal":null},"K":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.002925,"minimum":null,"nominal":null}},"family":"ep","familySubtype":null,"magneticCircuit":"open","name":"EP 17","type":"standard"},"type":"half set"},{"coordinates":[0.0,0.0,0.0],"dimensions":null,"insulationMaterial":null,"machining":[{"coordinates":[0.0,-0.0001035,0.0],"length":0.000207}],"material":"3C91","rotation":[0.0,0.0,0.0],"shape":{"aliases":[],"dimensions":{"A":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0184,"minimum":0.0176,"nominal":null},"B":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0085,"minimum":0.0083,"nominal":null},"C":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.01125,"minimum":0.01075,"nominal":null},"D":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0058,"minimum":0.0055,"nominal":null},"E":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.0124,"minimum":0.0116,"nominal":null},"F":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.00585,"minimum":0.0055,"nominal":null},"K":{"excludeMaximum":null,"excludeMinimum":null,"maximum":0.002925,"minimum":null,"nominal":null}},"family":"ep","familySubtype":null,"magneticCircuit":"open","name":"EP 17","type":"standard"},"type":"half set"}],"manufacturerInfo":{"cost":null,"datasheetUrl":"https://ferroxcube.com/upload/media/product/file/Pr_ds/EP17.pdf","family":null,"name":"Ferroxcube","orderCode":null,"reference":"EP17-3C91-A100","status":"production"},"name":"EP 17 - 3C91 - Gapped 0.414 mm","processedDescription":{"columns":[{"area":2.6e-05,"coordinates":[0.0,0.0,0.0],"depth":0.005675,"height":0.0113,"minimumDepth":null,"minimumWidth":null,"shape":"round","type":"central","width":0.005675},{"area":0.000163,"coordinates":[0.0,0.0,-0.007038],"depth":0.002075,"height":0.0113,"minimumDepth":null,"minimumWidth":0.003001,"shape":"irregular","type":"lateral","width":0.078555}],"depth":0.011,"effectiveParameters":{"effectiveArea":3.4461049794381896e-05,"effectiveLength":0.02849715978271948,"effectiveVolume":9.820420422707531e-07,"minimumArea":2.5790801226066944e-05},"height":0.016800000000000002,"width":0.018000000000000002,"windingWindows":[{"angle":null,"area":3.573625e-05,"coordinates":[0.0028374999999999997,0.0],"height":0.0113,"radialHeight":null,"sectionsOrientation":null,"shape":null,"width":0.0031625000000000004}]}},"distributorsInfo":null,"manufacturerInfo":null,"rotation":null},"outputs":[]})"));

    auto bobbin = mas.get_magnetic().get_coil().get_bobbin();
    if (std::holds_alternative<std::string>(bobbin)) {
        auto bobbinString = std::get<std::string>(bobbin);
        if (bobbinString == "Dummy") {
            mas.get_mutable_magnetic().get_mutable_coil().set_bobbin(OpenMagnetics::Bobbin::create_quick_bobbin(mas.get_mutable_magnetic().get_mutable_core()));
        }
    }
    mas.get_mutable_inputs().process();

    auto sections = CoilAdviser().get_advised_sections(mas, pattern, repetitions);
}

TEST_CASE("Test_CoilAdviser_Web_3", "[adviser][coil-adviser][bug][smoke-test]") {
    auto json_path_2049 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_coiladviser_web_3_2049.json");
    std::ifstream json_file_2049(json_path_2049);
    OpenMagnetics::Mas mas(json::parse(json_file_2049));


    for (size_t windingIndex = 0; windingIndex < mas.get_magnetic().get_coil().get_functional_description().size(); ++windingIndex) {
        mas.get_mutable_magnetic().get_mutable_coil().get_mutable_functional_description()[windingIndex].set_wire("Dummy");
    }
    mas.get_mutable_magnetic().get_mutable_coil().set_turns_description(std::nullopt);
    mas.get_mutable_magnetic().get_mutable_coil().set_layers_description(std::nullopt);
    mas.get_mutable_magnetic().get_mutable_coil().set_sections_description(std::nullopt);

    CoilAdviser coilAdviser;
    auto masMagneticsWithCoil = coilAdviser.get_advised_coil(mas, 1);
}

TEST_CASE("Test_CoilAdviser_Web_4", "[adviser][coil-adviser][bug][smoke-test]") {
    auto json_path_2064 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_coiladviser_web_4_2064.json");
    std::ifstream json_file_2064(json_path_2064);
    OpenMagnetics::Mas mas(json::parse(json_file_2064));


    for (size_t windingIndex = 0; windingIndex < mas.get_magnetic().get_coil().get_functional_description().size(); ++windingIndex) {
        mas.get_mutable_magnetic().get_mutable_coil().get_mutable_functional_description()[windingIndex].set_wire("Dummy");
    }
    mas.get_mutable_magnetic().get_mutable_coil().set_turns_description(std::nullopt);
    mas.get_mutable_magnetic().get_mutable_coil().set_layers_description(std::nullopt);
    mas.get_mutable_magnetic().get_mutable_coil().set_sections_description(std::nullopt);
    mas.get_mutable_magnetic().get_mutable_coil().set_groups_description(std::nullopt);

    CoilAdviser coilAdviser;
    auto mases = coilAdviser.get_advised_coil(mas, 1);
    REQUIRE(mases.size() > 0);

    auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
    auto outFile = outputFilePath;
    std::string filename = "Test_CoilAdviser_Web_4.svg";
    outFile.append(filename);
    Painter painter(outFile);

    painter.paint_core(mases[0].get_mutable_magnetic());
    painter.paint_bobbin(mases[0].get_mutable_magnetic());
    painter.paint_coil_turns(mases[0].get_mutable_magnetic());
    // painter.paint_coil_sections(mases[0].get_mutable_magnetic());
    painter.export_svg();
}

TEST_CASE("Test_CoilAdviser_Web_5", "[adviser][coil-adviser][bug][smoke-test]") {
    auto json_path_2093 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_coiladviser_web_5_2093.json");
    std::ifstream json_file_2093(json_path_2093);
    OpenMagnetics::Mas mas(json::parse(json_file_2093));


    for (size_t windingIndex = 0; windingIndex < mas.get_magnetic().get_coil().get_functional_description().size(); ++windingIndex) {
        mas.get_mutable_magnetic().get_mutable_coil().get_mutable_functional_description()[windingIndex].set_wire("Dummy");
    }
    mas.get_mutable_magnetic().get_mutable_coil().set_turns_description(std::nullopt);
    mas.get_mutable_magnetic().get_mutable_coil().set_layers_description(std::nullopt);
    mas.get_mutable_magnetic().get_mutable_coil().set_sections_description(std::nullopt);
    mas.get_mutable_magnetic().get_mutable_coil().set_groups_description(std::nullopt);

    CoilAdviser coilAdviser;
    auto mases = coilAdviser.get_advised_coil(mas, 1);
    REQUIRE(mases.size() > 0);

    auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
    auto outFile = outputFilePath;
    std::string filename = "Test_CoilAdviser_Web_5.svg";
    outFile.append(filename);
    Painter painter(outFile);

    painter.paint_core(mases[0].get_mutable_magnetic());
    painter.paint_bobbin(mases[0].get_mutable_magnetic());
    painter.paint_coil_turns(mases[0].get_mutable_magnetic());
    // painter.paint_coil_sections(mases[0].get_mutable_magnetic());
    painter.export_svg();
}

TEST_CASE("Test_CoilAdviser_Web_6", "[adviser][coil-adviser][bug][smoke-test]") {
    auto json_path_2122 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_coiladviser_web_6_2122.json");
    std::ifstream json_file_2122(json_path_2122);
    OpenMagnetics::Mas mas(json::parse(json_file_2122));


    for (size_t windingIndex = 0; windingIndex < mas.get_magnetic().get_coil().get_functional_description().size(); ++windingIndex) {
        mas.get_mutable_magnetic().get_mutable_coil().get_mutable_functional_description()[windingIndex].set_wire("Dummy");
    }
    mas.get_mutable_magnetic().get_mutable_coil().set_turns_description(std::nullopt);
    mas.get_mutable_magnetic().get_mutable_coil().set_layers_description(std::nullopt);
    mas.get_mutable_magnetic().get_mutable_coil().set_sections_description(std::nullopt);
    mas.get_mutable_magnetic().get_mutable_coil().set_groups_description(std::nullopt);

    CoilAdviser coilAdviser;
    auto mases = coilAdviser.get_advised_coil(mas, 1);
    REQUIRE(mases.size() > 0);
}

}  // namespace
