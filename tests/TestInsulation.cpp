#include "Insulation.h"
#include <UnitTest++.h>
#include "TestingUtils.h"


SUITE(Insulation) {
    TEST(IEC_60664) {
        auto standard = OpenMagnetics::InsulationIEC60664Model();
        OpenMagnetics::DimensionWithTolerance altitude;
        altitude.set_maximum(2000);
        auto cti = OpenMagnetics::Cti::GROUP_II;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        OpenMagnetics::DimensionWithTolerance mainSupplyVoltage;
        mainSupplyVoltage.set_nominal(400);
        auto overvoltageCategory = OpenMagnetics::OvervoltageCategory::OVC_II;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P2;
        auto standards = std::vector<OpenMagnetics::InsulationStandards>{OpenMagnetics::InsulationStandards::IEC_606641};
        double maximumVoltageRms = 666;
        double maximumVoltagePeak = 800;
        double frequency = 30000;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND); auto solidInsulation = standard.calculate_solid_insulation(inputs);
        std::cout << "solidInsulation: " << solidInsulation << std::endl;
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        std::cout << "creepageDistance: " << creepageDistance << std::endl;
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
        frequency = 800000;
        altitude.set_maximum(2000);
        mainSupplyVoltage.set_nominal(400);
        auto cti = OpenMagnetics::Cti::GROUP_I;
        auto insulationType = OpenMagnetics::InsulationType::BASIC;
        auto pollutionDegree = OpenMagnetics::PollutionDegree::P1;
        OpenMagnetics::InputsWrapper inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, cti, insulationType, mainSupplyVoltage, overvoltageCategory, pollutionDegree, standards, maximumVoltageRms, maximumVoltagePeak, frequency, OpenMagnetics::WiringTechnology::WOUND);
        auto creepageDistance = standard.calculate_creepage_distance(inputs);
        CHECK(0.0024 < creepageDistance);
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
        CHECK(0.0048 < creepageDistance);
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
        CHECK(0.004 < creepageDistance);
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
        CHECK(0.008 < creepageDistance);
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
        CHECK(0.01 < creepageDistance);
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
        CHECK(0.02 < creepageDistance);
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
        CHECK(0.0024 < creepageDistance);
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
        CHECK(0.0048 < creepageDistance);
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
        CHECK(0.0056 < creepageDistance);
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
        CHECK(0.0112 < creepageDistance);
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
        CHECK(0.011 < creepageDistance);
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
        CHECK(0.022 < creepageDistance);
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
        CHECK(0.0024 < creepageDistance);
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
        CHECK(0.0048 < creepageDistance);
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
        CHECK(0.008 < creepageDistance);
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
        CHECK(0.016 < creepageDistance);
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
        CHECK(0.0125 < creepageDistance);
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
        CHECK(0.025 < creepageDistance);
    }

    // testear rangos de frequencia creepage
    // poluciones creepage
    // reinforced y double
    // altitudes

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
        CHECK_EQUAL(0.00006, clearance);
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
        CHECK_EQUAL(0.00015, clearance);
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
        CHECK_EQUAL(0.0002, clearance);
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
        CHECK_EQUAL(0.001, clearance);
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
        CHECK_EQUAL(0.001, clearance);
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
        CHECK_EQUAL(0.001, clearance);
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
        CHECK_EQUAL(0.002, clearance);
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
        CHECK_EQUAL(0.002, clearance);
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
        CHECK_EQUAL(0.002, clearance);
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