#include "MAS.hpp"
#include "constructive_models/Wire.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/Utils.h"

#include "TestingUtils.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace MAS;
using namespace OpenMagnetics;

namespace {
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

    TEST_CASE("Test_Round_Wire_20C", "[physical-model][ohmic-losses][smoke-test]") {
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
        REQUIRE_THAT(expectedDcResistance, Catch::Matchers::WithinAbs(dcResistance, expectedDcResistance * maximumError));
    }

    TEST_CASE("Test_Round_Wire_200C", "[physical-model][ohmic-losses][smoke-test]") {
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
        REQUIRE_THAT(expectedDcResistance, Catch::Matchers::WithinAbs(dcResistance, expectedDcResistance * maximumError));
    }

    TEST_CASE("Test_Litz_Wire_Small", "[physical-model][ohmic-losses][smoke-test]") {
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
        REQUIRE_THAT(expectedDcResistance, Catch::Matchers::WithinAbs(dcResistance, expectedDcResistance * maximumError));
    }

    TEST_CASE("Test_Litz_Wire_Large", "[physical-model][ohmic-losses][smoke-test]") {
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
        REQUIRE_THAT(expectedDcResistance, Catch::Matchers::WithinAbs(dcResistance, expectedDcResistance * maximumError));
    }

    TEST_CASE("Test_Foil_Wire_20C", "[physical-model][ohmic-losses][smoke-test]") {
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
        REQUIRE_THAT(expectedDcResistance, Catch::Matchers::WithinAbs(dcResistance, expectedDcResistance * maximumError));
    }

    TEST_CASE("Test_Rectangular_Wire_20C", "[physical-model][ohmic-losses][smoke-test]") {
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
        REQUIRE_THAT(expectedDcResistance, Catch::Matchers::WithinAbs(dcResistance, expectedDcResistance * maximumError));
    }

    TEST_CASE("Test_Winding_Ohmic_Losses_One_Turn", "[physical-model][ohmic-losses][smoke-test]") {
        double temperature = 20;
        std::vector<int64_t> numberTurns = {1};
        std::vector<int64_t> numberParallels = {1};

        auto operatingPoint = get_operating_point_with_dc_current({1});
        auto winding = get_coil(numberTurns, numberParallels);

        auto windingOhmicLosses = WindingOhmicLosses();
        double ohmicLosses = windingOhmicLosses.calculate_ohmic_losses(winding, operatingPoint, temperature).get_winding_losses();
        double expectedOhmicLosses = 3.1e-3;
        REQUIRE_THAT(expectedOhmicLosses, Catch::Matchers::WithinAbs(ohmicLosses, expectedOhmicLosses * maximumError));
    }

    TEST_CASE("Test_Winding_Ohmic_Losses_Two_Turns", "[physical-model][ohmic-losses][smoke-test]") {
        double temperature = 20;
        std::vector<int64_t> numberTurns = {2};
        std::vector<int64_t> numberParallels = {1};

        auto operatingPoint = get_operating_point_with_dc_current({1});
        auto winding = get_coil(numberTurns, numberParallels);

        auto windingOhmicLosses = WindingOhmicLosses();
        double ohmicLosses = windingOhmicLosses.calculate_ohmic_losses(winding, operatingPoint, temperature).get_winding_losses();
        double expectedOhmicLosses = 6.2e-3;
        REQUIRE_THAT(expectedOhmicLosses, Catch::Matchers::WithinAbs(ohmicLosses, expectedOhmicLosses * maximumError));
    }

    TEST_CASE("Test_Winding_Ohmic_Losses_Two_Turns_Two_Parallels", "[physical-model][ohmic-losses][smoke-test]") {
        double temperature = 20;
        std::vector<int64_t> numberTurns = {2};
        std::vector<int64_t> numberParallels = {2};

        auto operatingPoint = get_operating_point_with_dc_current({1});
        auto winding = get_coil(numberTurns, numberParallels);

        auto windingOhmicLosses = WindingOhmicLosses();
        double ohmicLosses = windingOhmicLosses.calculate_ohmic_losses(winding, operatingPoint, temperature).get_winding_losses();
        double expectedOhmicLosses = 3.1e-3;
        REQUIRE_THAT(expectedOhmicLosses, Catch::Matchers::WithinAbs(ohmicLosses, expectedOhmicLosses * maximumError));
    }

    TEST_CASE("Test_Winding_Ohmic_Losses_One_Turn_Double_Current", "[physical-model][ohmic-losses][smoke-test]") {
        double temperature = 20;
        std::vector<int64_t> numberTurns = {1};
        std::vector<int64_t> numberParallels = {1};

        auto operatingPoint = get_operating_point_with_dc_current({2});
        auto winding = get_coil(numberTurns, numberParallels);

        auto windingOhmicLosses = WindingOhmicLosses();
        double ohmicLosses = windingOhmicLosses.calculate_ohmic_losses(winding, operatingPoint, temperature).get_winding_losses();
        double expectedOhmicLosses = 12.4e-3;
        REQUIRE_THAT(expectedOhmicLosses, Catch::Matchers::WithinAbs(ohmicLosses, expectedOhmicLosses * maximumError));
    }

    TEST_CASE("Test_Winding_Ohmic_Losses_Two_Windings", "[physical-model][ohmic-losses][smoke-test]") {
        double temperature = 20;
        std::vector<int64_t> numberTurns = {1, 2};
        std::vector<int64_t> numberParallels = {1, 2};

        auto operatingPoint = get_operating_point_with_dc_current({1, 1});
        auto winding = get_coil(numberTurns, numberParallels);

        auto windingOhmicLosses = WindingOhmicLosses();
        double ohmicLosses = windingOhmicLosses.calculate_ohmic_losses(winding, operatingPoint, temperature).get_winding_losses();
        double expectedOhmicLosses = 6.55e-3;
        REQUIRE_THAT(expectedOhmicLosses, Catch::Matchers::WithinAbs(ohmicLosses, expectedOhmicLosses * maximumError));
    }

    TEST_CASE("Test_Winding_Ohmic_Losses_Two_Windings_Double_Turns", "[physical-model][ohmic-losses][smoke-test]") {
        double temperature = 20;
        std::vector<int64_t> numberTurns = {2, 4};
        std::vector<int64_t> numberParallels = {1, 2};

        auto operatingPoint = get_operating_point_with_dc_current({1, 1});
        auto winding = get_coil(numberTurns, numberParallels);

        auto windingOhmicLosses = WindingOhmicLosses();
        double ohmicLosses = windingOhmicLosses.calculate_ohmic_losses(winding, operatingPoint, temperature).get_winding_losses();
        double expectedOhmicLosses = 2 * 6.55e-3;
        REQUIRE_THAT(expectedOhmicLosses, Catch::Matchers::WithinAbs(ohmicLosses, expectedOhmicLosses * maximumError));
    }

    TEST_CASE("Test_Winding_Ohmic_Losses_Two_Windings_High_Temp", "[physical-model][ohmic-losses][smoke-test]") {
        double temperature = 120;
        std::vector<int64_t> numberTurns = {1, 2};
        std::vector<int64_t> numberParallels = {1, 2};

        auto operatingPoint = get_operating_point_with_dc_current({1, 1});
        auto winding = get_coil(numberTurns, numberParallels);

        auto windingOhmicLosses = WindingOhmicLosses();
        double ohmicLosses = windingOhmicLosses.calculate_ohmic_losses(winding, operatingPoint, temperature).get_winding_losses();
        double expectedOhmicLosses = 9.2e-3;
        REQUIRE_THAT(expectedOhmicLosses, Catch::Matchers::WithinAbs(ohmicLosses, expectedOhmicLosses * maximumError));
    }


}  // namespace
