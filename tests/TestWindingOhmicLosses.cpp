#include "MAS.hpp"
#include "constructive_models/Wire.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/Utils.h"

#include "TestingUtils.h"
#include <UnitTest++.h>

using namespace MAS;
using namespace OpenMagnetics;

SUITE(WindingOhmicLosses) {
    double maximumError = 0.05;

    OperatingPoint get_operating_point_with_dc_current(std::vector<double> dcCurrents) {
        OperatingPoint operatingPoint;
        OperatingPointExcitation operatingPointExcitation;
        SignalDescriptor current;
        Processed processed;
        std::vector<OperatingPointExcitation> excitations;
        for (auto& dcCurrent : dcCurrents) {
            processed.set_rms(dcCurrent);
            current.set_processed(processed);
            operatingPointExcitation.set_current(current);
            excitations.push_back(operatingPointExcitation);
        }
        operatingPoint.set_excitations_per_winding(excitations);

        return operatingPoint;
    }

    OpenMagnetics::Coil get_coil(std::vector<int64_t> numberTurns, std::vector<int64_t> numberParallels) {
        int64_t numberPhysicalTurns = std::numeric_limits<int64_t>::max();
        for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex)
        {
            numberPhysicalTurns = std::min(numberPhysicalTurns, numberTurns[windingIndex] * numberParallels[windingIndex]);
        }
        double bobbinHeight = 0.01;
        double bobbinWidth = 0.01;
        std::vector<double> bobbinCenterCoodinates = {0.01, 0, 0};
        uint8_t interleavingLevel = 1;
        interleavingLevel = std::min(uint8_t(numberPhysicalTurns), interleavingLevel);
        auto windingOrientation = WindingOrientation::OVERLAPPING;

        auto winding = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, bobbinHeight, bobbinWidth, bobbinCenterCoodinates, interleavingLevel, windingOrientation);
        return winding;
    }

    TEST(Test_Round_Wire_20C) {
        double temperature = 20;
        Turn turn;
        turn.set_length(1);

        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_diameter(0.00032114);
        wire.set_material("copper");
        wire.set_type(WireType::ROUND);

        auto windingOhmicLosses = WindingOhmicLosses();
        double dcResistance = windingOhmicLosses.calculate_dc_resistance(turn, wire, temperature);
        double expectedDcResistance = 211.1e-3;
        CHECK_CLOSE(expectedDcResistance, dcResistance, expectedDcResistance * maximumError);
    }

    TEST(Test_Round_Wire_200C) {
        double temperature = 200;
        Turn turn;
        turn.set_length(1);

        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_diameter(0.00032114);
        wire.set_material("copper");
        wire.set_type(WireType::ROUND);

        auto windingOhmicLosses = WindingOhmicLosses();
        double dcResistance = windingOhmicLosses.calculate_dc_resistance(turn, wire, temperature);
        double expectedDcResistance = 357e-3;
        CHECK_CLOSE(expectedDcResistance, dcResistance, expectedDcResistance * maximumError);
    }

    TEST(Test_Litz_Wire_Small) {
        double temperature = 20;
        Turn turn;
        turn.set_length(1);

        OpenMagnetics::Wire wire;
        WireRound strand;
        strand.set_type(WireType::ROUND);
        strand.set_material("copper");
        DimensionWithTolerance conductingDiameter;
        conductingDiameter.set_nominal(0.000040);
        strand.set_conducting_diameter(conductingDiameter);
        wire.set_strand(strand);
        wire.set_number_conductors(30);
        wire.set_type(WireType::LITZ);

        auto windingOhmicLosses = WindingOhmicLosses();
        double dcResistance = windingOhmicLosses.calculate_dc_resistance(turn, wire, temperature);
        double expectedDcResistance = 0.4625;
        CHECK_CLOSE(expectedDcResistance, dcResistance, expectedDcResistance * maximumError);
    }

    TEST(Test_Litz_Wire_Large) {
        double temperature = 20;
        Turn turn;
        turn.set_length(1);

        OpenMagnetics::Wire wire;
        WireRound strand;
        strand.set_type(WireType::ROUND);
        strand.set_material("copper");
        DimensionWithTolerance conductingDiameter;
        conductingDiameter.set_nominal(0.00012);
        strand.set_conducting_diameter(conductingDiameter);
        wire.set_strand(strand);
        wire.set_number_conductors(600);
        wire.set_type(WireType::LITZ);

        auto windingOhmicLosses = WindingOhmicLosses();
        double dcResistance = windingOhmicLosses.calculate_dc_resistance(turn, wire, temperature);
        double expectedDcResistance = 0.0025;
        CHECK_CLOSE(expectedDcResistance, dcResistance, expectedDcResistance * maximumError);
    }

    TEST(Test_Foil_Wire_20C) {
        double temperature = 20;
        Turn turn;
        turn.set_length(1);

        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_width(0.0005);
        wire.set_nominal_value_conducting_height(0.01);
        wire.set_material("copper");
        wire.set_type(WireType::FOIL);

        auto windingOhmicLosses = WindingOhmicLosses();
        double dcResistance = windingOhmicLosses.calculate_dc_resistance(turn, wire, temperature);
        double expectedDcResistance = 3.3e-3;
        CHECK_CLOSE(expectedDcResistance, dcResistance, expectedDcResistance * maximumError);
    }

    TEST(Test_Rectangular_Wire_20C) {
        double temperature = 20;
        Turn turn;
        turn.set_length(1);

        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_width(0.005);
        wire.set_nominal_value_conducting_height(0.001);
        wire.set_material("copper");
        wire.set_type(WireType::FOIL);

        auto windingOhmicLosses = WindingOhmicLosses();
        double dcResistance = windingOhmicLosses.calculate_dc_resistance(turn, wire, temperature);
        double expectedDcResistance = 3.3e-3;
        CHECK_CLOSE(expectedDcResistance, dcResistance, expectedDcResistance * maximumError);
    }

    TEST(Test_Winding_Ohmic_Losses_One_Turn) {
        double temperature = 20;
        std::vector<int64_t> numberTurns = {1};
        std::vector<int64_t> numberParallels = {1};

        auto operatingPoint = get_operating_point_with_dc_current({1});
        auto winding = get_coil(numberTurns, numberParallels);

        auto windingOhmicLosses = WindingOhmicLosses();
        double ohmicLosses = windingOhmicLosses.calculate_ohmic_losses(winding, operatingPoint, temperature).get_winding_losses();
        double expectedOhmicLosses = 3.1e-3;
        CHECK_CLOSE(expectedOhmicLosses, ohmicLosses, expectedOhmicLosses * maximumError);
    }

    TEST(Test_Winding_Ohmic_Losses_Two_Turns) {
        double temperature = 20;
        std::vector<int64_t> numberTurns = {2};
        std::vector<int64_t> numberParallels = {1};

        auto operatingPoint = get_operating_point_with_dc_current({1});
        auto winding = get_coil(numberTurns, numberParallels);

        auto windingOhmicLosses = WindingOhmicLosses();
        double ohmicLosses = windingOhmicLosses.calculate_ohmic_losses(winding, operatingPoint, temperature).get_winding_losses();
        double expectedOhmicLosses = 6.2e-3;
        CHECK_CLOSE(expectedOhmicLosses, ohmicLosses, expectedOhmicLosses * maximumError);
    }

    TEST(Test_Winding_Ohmic_Losses_Two_Turns_Two_Parallels) {
        double temperature = 20;
        std::vector<int64_t> numberTurns = {2};
        std::vector<int64_t> numberParallels = {2};

        auto operatingPoint = get_operating_point_with_dc_current({1});
        auto winding = get_coil(numberTurns, numberParallels);

        auto windingOhmicLosses = WindingOhmicLosses();
        double ohmicLosses = windingOhmicLosses.calculate_ohmic_losses(winding, operatingPoint, temperature).get_winding_losses();
        double expectedOhmicLosses = 3.1e-3;
        CHECK_CLOSE(expectedOhmicLosses, ohmicLosses, expectedOhmicLosses * maximumError);
    }

    TEST(Test_Winding_Ohmic_Losses_One_Turn_Double_Current) {
        double temperature = 20;
        std::vector<int64_t> numberTurns = {1};
        std::vector<int64_t> numberParallels = {1};

        auto operatingPoint = get_operating_point_with_dc_current({2});
        auto winding = get_coil(numberTurns, numberParallels);

        auto windingOhmicLosses = WindingOhmicLosses();
        double ohmicLosses = windingOhmicLosses.calculate_ohmic_losses(winding, operatingPoint, temperature).get_winding_losses();
        double expectedOhmicLosses = 12.4e-3;
        CHECK_CLOSE(expectedOhmicLosses, ohmicLosses, expectedOhmicLosses * maximumError);
    }

    TEST(Test_Winding_Ohmic_Losses_Two_Windings) {
        double temperature = 20;
        std::vector<int64_t> numberTurns = {1, 2};
        std::vector<int64_t> numberParallels = {1, 2};

        auto operatingPoint = get_operating_point_with_dc_current({1, 1});
        auto winding = get_coil(numberTurns, numberParallels);

        auto windingOhmicLosses = WindingOhmicLosses();
        double ohmicLosses = windingOhmicLosses.calculate_ohmic_losses(winding, operatingPoint, temperature).get_winding_losses();
        double expectedOhmicLosses = 6.55e-3;
        CHECK_CLOSE(expectedOhmicLosses, ohmicLosses, expectedOhmicLosses * maximumError);
    }

    TEST(Test_Winding_Ohmic_Losses_Two_Windings_Double_Turns) {
        double temperature = 20;
        std::vector<int64_t> numberTurns = {2, 4};
        std::vector<int64_t> numberParallels = {1, 2};

        auto operatingPoint = get_operating_point_with_dc_current({1, 1});
        auto winding = get_coil(numberTurns, numberParallels);

        auto windingOhmicLosses = WindingOhmicLosses();
        double ohmicLosses = windingOhmicLosses.calculate_ohmic_losses(winding, operatingPoint, temperature).get_winding_losses();
        double expectedOhmicLosses = 2 * 6.55e-3;
        CHECK_CLOSE(expectedOhmicLosses, ohmicLosses, expectedOhmicLosses * maximumError);
    }

    TEST(Test_Winding_Ohmic_Losses_Two_Windings_High_Temp) {
        double temperature = 120;
        std::vector<int64_t> numberTurns = {1, 2};
        std::vector<int64_t> numberParallels = {1, 2};

        auto operatingPoint = get_operating_point_with_dc_current({1, 1});
        auto winding = get_coil(numberTurns, numberParallels);

        auto windingOhmicLosses = WindingOhmicLosses();
        double ohmicLosses = windingOhmicLosses.calculate_ohmic_losses(winding, operatingPoint, temperature).get_winding_losses();
        double expectedOhmicLosses = 9.2e-3;
        CHECK_CLOSE(expectedOhmicLosses, ohmicLosses, expectedOhmicLosses * maximumError);
    }

}