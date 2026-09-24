#include "RandomUtils.h"
#include "constructive_models/NumberTurns.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;

namespace { 
    TEST_CASE("Number_Turns_Inductor", "[constructive-model][number-turns][smoke-test]") {
        DesignRequirements designRequirements;
        designRequirements.set_turns_ratios(std::vector<DimensionWithTolerance>{});
        uint64_t initialPrimaryNumberTurns = 42;

        NumberTurns numberTurns(initialPrimaryNumberTurns, designRequirements);
        std::vector<uint64_t> numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        REQUIRE(numberTurnsCombination[0] == initialPrimaryNumberTurns);
        numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        REQUIRE(numberTurnsCombination[0] == initialPrimaryNumberTurns + 1);
    }

    TEST_CASE("Number_Turns_Two_Windings_Turns_Ratio_1", "[constructive-model][number-turns][smoke-test]") {
        DesignRequirements designRequirements;
        DimensionWithTolerance turnsRatio;
        double turnsRatioValue = 1;
        turnsRatio.set_nominal(turnsRatioValue);
        designRequirements.set_turns_ratios(std::vector<DimensionWithTolerance>{turnsRatio});
        uint64_t initialPrimaryNumberTurns = 42;

        NumberTurns numberTurns(initialPrimaryNumberTurns, designRequirements);
        std::vector<uint64_t> numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        REQUIRE(numberTurnsCombination[0] == initialPrimaryNumberTurns);
        REQUIRE(numberTurnsCombination[1] == initialPrimaryNumberTurns * 1);
        numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        REQUIRE(numberTurnsCombination[0] == initialPrimaryNumberTurns + 1);
        REQUIRE(numberTurnsCombination[1] == (initialPrimaryNumberTurns + 1) * 1);
    }

    // A bound is inclusive unless MAS marks it excluded. The one-sided branches of
    // check_requirement were strict, so a value exactly at a maximum-only requirement
    // failed while the same value passed a two-sided one.
    TEST_CASE("Check_Requirement_One_Sided_Bounds_Are_Inclusive_Unless_Excluded", "[support][check-requirement]") {
        DimensionWithTolerance ceiling;
        ceiling.set_maximum(0.5);
        CHECK(check_requirement(ceiling, 0.5));
        CHECK(check_requirement(ceiling, 0.4));
        CHECK_FALSE(check_requirement(ceiling, 0.51));
        ceiling.set_exclude_maximum(true);
        CHECK_FALSE(check_requirement(ceiling, 0.5));
        CHECK(check_requirement(ceiling, 0.49));

        DimensionWithTolerance floor;
        floor.set_minimum(70e-6);
        CHECK(check_requirement(floor, 70e-6));
        CHECK_FALSE(check_requirement(floor, 69e-6));
        floor.set_exclude_minimum(true);
        CHECK_FALSE(check_requirement(floor, 70e-6));

        DimensionWithTolerance band;
        band.set_minimum(1.0);
        band.set_maximum(2.0);
        band.set_exclude_maximum(true);
        CHECK(check_requirement(band, 1.0));
        CHECK_FALSE(check_requirement(band, 2.0));
    }

    // The Weinberg transformer: a centre-tapped primary (1:1) and two secondaries whose
    // ratio is a CEILING (<= 0.5, i.e. at least twice the primary turns). NumberTurns aims
    // each ceiling exactly at its maximum; with the strict one-sided check that exact
    // value never passed and every call threw "NumberTurns did not converge", so the
    // core adviser culled every candidate and advised nothing.
    TEST_CASE("Number_Turns_Converges_On_Maximum_Only_Turns_Ratios", "[constructive-model][number-turns]") {
        DesignRequirements designRequirements;
        DimensionWithTolerance primaryHalf;
        primaryHalf.set_nominal(1.0);
        DimensionWithTolerance secondaryCeiling;
        secondaryCeiling.set_maximum(0.5);
        designRequirements.set_turns_ratios(std::vector<DimensionWithTolerance>{primaryHalf, secondaryCeiling, secondaryCeiling});

        NumberTurns numberTurns(7, designRequirements);
        auto combination = numberTurns.get_next_number_turns_combination();
        REQUIRE(combination.size() == 4);
        CHECK(combination[1] == combination[0]);
        CHECK(check_requirement(secondaryCeiling, double(combination[0]) / combination[2]));
        CHECK(check_requirement(secondaryCeiling, double(combination[0]) / combination[3]));
    }

    TEST_CASE("Number_Turns_Two_Windings_Turns_Ratio_8", "[constructive-model][number-turns][smoke-test]") {
        DesignRequirements designRequirements;
        DimensionWithTolerance turnsRatio;
        double turnsRatioValue = 8;
        turnsRatio.set_nominal(turnsRatioValue);
        turnsRatio.set_minimum(turnsRatioValue * 0.8);
        turnsRatio.set_maximum(turnsRatioValue * 1.2);
        designRequirements.set_turns_ratios(std::vector<DimensionWithTolerance>{turnsRatio});
        uint64_t initialPrimaryNumberTurns = 42;

        NumberTurns numberTurns(initialPrimaryNumberTurns, designRequirements);
        std::vector<uint64_t> numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        REQUIRE(numberTurnsCombination[0] == initialPrimaryNumberTurns);
        REQUIRE(check_requirement(turnsRatio, double(numberTurnsCombination[0]) / numberTurnsCombination[1]));

        numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        REQUIRE(numberTurnsCombination[0] == initialPrimaryNumberTurns + 1);
        REQUIRE(check_requirement(turnsRatio, double(numberTurnsCombination[0]) / numberTurnsCombination[1]));
    }

    TEST_CASE("Number_Turns_Two_Windings_Turns_Ratio_0_001", "[constructive-model][number-turns][smoke-test]") {
        DesignRequirements designRequirements;
        DimensionWithTolerance turnsRatio;
        double turnsRatioValue = 0.001;
        turnsRatio.set_nominal(turnsRatioValue);
        turnsRatio.set_minimum(turnsRatioValue * 0.8);
        turnsRatio.set_maximum(turnsRatioValue * 1.2);
        designRequirements.set_turns_ratios(std::vector<DimensionWithTolerance>{turnsRatio});
        uint64_t initialPrimaryNumberTurns = 42;

        NumberTurns numberTurns(initialPrimaryNumberTurns, designRequirements);
        std::vector<uint64_t> numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        REQUIRE(numberTurnsCombination[0] == initialPrimaryNumberTurns);
        REQUIRE(check_requirement(turnsRatio, double(numberTurnsCombination[0]) / numberTurnsCombination[1]));

        numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        REQUIRE(numberTurnsCombination[0] == initialPrimaryNumberTurns + 1);
        REQUIRE(check_requirement(turnsRatio, double(numberTurnsCombination[0]) / numberTurnsCombination[1]));
    }

    TEST_CASE("Number_Turns_Two_Windings_Turns_Ratio_Random", "[constructive-model][number-turns]") {
        for (size_t i = 0; i < 1000; ++i)
        {
            DesignRequirements designRequirements;
            DimensionWithTolerance turnsRatio;
            double turnsRatioValue =  ((double) OpenMagnetics::TestUtils::randomInt(0, RAND_MAX) / RAND_MAX) * (100 - 0.0001) + 0.0001;
            if (OpenMagnetics::TestUtils::randomInt(0, 2 - 1) == 0) {
                turnsRatioValue = 1 / turnsRatioValue;
            }
            turnsRatio.set_nominal(turnsRatioValue);
            turnsRatio.set_minimum(turnsRatioValue * 0.95);
            turnsRatio.set_maximum(turnsRatioValue * 1.05);
            designRequirements.set_turns_ratios(std::vector<DimensionWithTolerance>{turnsRatio});
            uint64_t initialPrimaryNumberTurns = OpenMagnetics::TestUtils::randomSize(1, 100 + 1 - 1);
            
            NumberTurns numberTurns(initialPrimaryNumberTurns, designRequirements);
            std::vector<uint64_t> numberTurnsCombination = numberTurns.get_next_number_turns_combination();
            REQUIRE(numberTurnsCombination[0] >= initialPrimaryNumberTurns);
            REQUIRE(check_requirement(turnsRatio, double(numberTurnsCombination[0]) / numberTurnsCombination[1]));

            numberTurnsCombination = numberTurns.get_next_number_turns_combination();
            REQUIRE(check_requirement(turnsRatio, double(numberTurnsCombination[0]) / numberTurnsCombination[1]));
        }
    }

    TEST_CASE("Number_Turns_Many_Windings_Turns_Ratio_Random", "[constructive-model][number-turns]") {
        for (size_t i = 0; i < 1000; ++i)
        {
            DesignRequirements designRequirements;
            std::vector<DimensionWithTolerance> turnsRatios;
            size_t numberSecondaryWindings = OpenMagnetics::TestUtils::randomInt(0, 10 - 1);
            for (size_t turnRatioIndex = 0; turnRatioIndex < numberSecondaryWindings; ++turnRatioIndex) {
                DimensionWithTolerance turnsRatio;
                double turnsRatioValue =  ((double) OpenMagnetics::TestUtils::randomInt(0, RAND_MAX) / RAND_MAX) * (100 - 0.0001) + 0.0001;
                if (OpenMagnetics::TestUtils::randomInt(0, 2 - 1) == 0) {
                    turnsRatioValue = 1 / turnsRatioValue;
                }
                turnsRatio.set_nominal(turnsRatioValue);
                turnsRatio.set_minimum(turnsRatioValue * 0.95);
                turnsRatio.set_maximum(turnsRatioValue * 1.05);
                turnsRatios.push_back(turnsRatio);

            }
            designRequirements.set_turns_ratios(turnsRatios);
            uint64_t initialPrimaryNumberTurns = OpenMagnetics::TestUtils::randomSize(1, 100 + 1 - 1);
            
            NumberTurns numberTurns(initialPrimaryNumberTurns, designRequirements);
            std::vector<uint64_t> numberTurnsCombination = numberTurns.get_next_number_turns_combination();
            REQUIRE(numberTurnsCombination[0] >= initialPrimaryNumberTurns);
            for (size_t turnRatioIndex = 0; turnRatioIndex < turnsRatios.size(); ++turnRatioIndex)
            {
                REQUIRE(check_requirement(turnsRatios[turnRatioIndex], double(numberTurnsCombination[0]) / numberTurnsCombination[turnRatioIndex + 1]));
            }

            numberTurnsCombination = numberTurns.get_next_number_turns_combination();
            for (size_t turnRatioIndex = 0; turnRatioIndex < turnsRatios.size(); ++turnRatioIndex)
            {
                REQUIRE(check_requirement(turnsRatios[turnRatioIndex], double(numberTurnsCombination[0]) / numberTurnsCombination[turnRatioIndex + 1]));
            }

        }
    }

    TEST_CASE("Number_Turns_Two_Windings_Turns_Ratio_Random_0", "[constructive-model][number-turns][smoke-test]") {
        DesignRequirements designRequirements;
        DimensionWithTolerance turnsRatio;
        double turnsRatioValue = 78;
        turnsRatio.set_nominal(turnsRatioValue);
        turnsRatio.set_minimum(turnsRatioValue * 0.8);
        turnsRatio.set_maximum(turnsRatioValue * 1.2);
        designRequirements.set_turns_ratios(std::vector<DimensionWithTolerance>{turnsRatio});
        uint64_t initialPrimaryNumberTurns = 40;
        
        NumberTurns numberTurns(initialPrimaryNumberTurns, designRequirements);
        std::vector<uint64_t> numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        REQUIRE(numberTurnsCombination[0] >= initialPrimaryNumberTurns);
        REQUIRE(check_requirement(turnsRatio, double(numberTurnsCombination[0]) / numberTurnsCombination[1]));

        numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        REQUIRE(check_requirement(turnsRatio, double(numberTurnsCombination[0]) / numberTurnsCombination[1]));
    }

    TEST_CASE("Number_Turns_Two_Windings_Turns_Ratio_Random_1", "[constructive-model][number-turns][smoke-test]") {
        DesignRequirements designRequirements;
        DimensionWithTolerance turnsRatio;
        double turnsRatioValue = 0.010101;
        turnsRatio.set_nominal(turnsRatioValue);
        turnsRatio.set_minimum(turnsRatioValue * 0.8);
        turnsRatio.set_maximum(turnsRatioValue * 1.2);
        designRequirements.set_turns_ratios(std::vector<DimensionWithTolerance>{turnsRatio});
        uint64_t initialPrimaryNumberTurns = 60;
        
        NumberTurns numberTurns(initialPrimaryNumberTurns, designRequirements);
        std::vector<uint64_t> numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        REQUIRE(numberTurnsCombination[0] >= initialPrimaryNumberTurns);
        REQUIRE(check_requirement(turnsRatio, double(numberTurnsCombination[0]) / numberTurnsCombination[1]));

        numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        REQUIRE(check_requirement(turnsRatio, double(numberTurnsCombination[0]) / numberTurnsCombination[1]));
    }

// End of SUITE

}  // namespace
