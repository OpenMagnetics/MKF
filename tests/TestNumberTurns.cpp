#include "constructive_models/NumberTurns.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;

namespace { 
    TEST_CASE("Number_Turns_Inductor", "[constructive-model][number-turns]") {
        DesignRequirements designRequirements;
        designRequirements.set_turns_ratios(std::vector<DimensionWithTolerance>{});
        uint64_t initialPrimaryNumberTurns = 42;

        NumberTurns numberTurns(initialPrimaryNumberTurns, designRequirements);
        std::vector<uint64_t> numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        REQUIRE(numberTurnsCombination[0] == initialPrimaryNumberTurns);
        numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        REQUIRE(numberTurnsCombination[0] == initialPrimaryNumberTurns + 1);
    }

    TEST_CASE("Number_Turns_Two_Windings_Turns_Ratio_1", "[constructive-model][number-turns]") {
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

    TEST_CASE("Number_Turns_Two_Windings_Turns_Ratio_8", "[constructive-model][number-turns]") {
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

    TEST_CASE("Number_Turns_Two_Windings_Turns_Ratio_0_001", "[constructive-model][number-turns]") {
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
        srand (time(NULL));
        for (size_t i = 0; i < 1000; ++i)
        {
            DesignRequirements designRequirements;
            DimensionWithTolerance turnsRatio;
            double turnsRatioValue =  ((double) std::rand() / RAND_MAX) * (100 - 0.0001) + 0.0001;
            if (std::rand() % 2 == 0) {
                turnsRatioValue = 1 / turnsRatioValue;
            }
            turnsRatio.set_nominal(turnsRatioValue);
            turnsRatio.set_minimum(turnsRatioValue * 0.95);
            turnsRatio.set_maximum(turnsRatioValue * 1.05);
            designRequirements.set_turns_ratios(std::vector<DimensionWithTolerance>{turnsRatio});
            uint64_t initialPrimaryNumberTurns = std::rand() % 100 + 1UL;
            
            NumberTurns numberTurns(initialPrimaryNumberTurns, designRequirements);
            std::vector<uint64_t> numberTurnsCombination = numberTurns.get_next_number_turns_combination();
            REQUIRE(numberTurnsCombination[0] >= initialPrimaryNumberTurns);
            REQUIRE(check_requirement(turnsRatio, double(numberTurnsCombination[0]) / numberTurnsCombination[1]));

            numberTurnsCombination = numberTurns.get_next_number_turns_combination();
            REQUIRE(check_requirement(turnsRatio, double(numberTurnsCombination[0]) / numberTurnsCombination[1]));
        }
    }

    TEST_CASE("Number_Turns_Many_Windings_Turns_Ratio_Random", "[constructive-model][number-turns]") {
        srand (time(NULL));
        for (size_t i = 0; i < 1000; ++i)
        {
            DesignRequirements designRequirements;
            std::vector<DimensionWithTolerance> turnsRatios;
            size_t numberSecondaryWindings = std::rand() % 10;
            for (size_t turnRatioIndex = 0; turnRatioIndex < numberSecondaryWindings; ++turnRatioIndex) {
                DimensionWithTolerance turnsRatio;
                double turnsRatioValue =  ((double) std::rand() / RAND_MAX) * (100 - 0.0001) + 0.0001;
                if (std::rand() % 2 == 0) {
                    turnsRatioValue = 1 / turnsRatioValue;
                }
                turnsRatio.set_nominal(turnsRatioValue);
                turnsRatio.set_minimum(turnsRatioValue * 0.95);
                turnsRatio.set_maximum(turnsRatioValue * 1.05);
                turnsRatios.push_back(turnsRatio);

            }
            designRequirements.set_turns_ratios(turnsRatios);
            uint64_t initialPrimaryNumberTurns = std::rand() % 100 + 1UL;
            
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

    TEST_CASE("Number_Turns_Two_Windings_Turns_Ratio_Random_0", "[constructive-model][number-turns]") {
        srand (time(NULL));
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

    TEST_CASE("Number_Turns_Two_Windings_Turns_Ratio_Random_1", "[constructive-model][number-turns]") {
        srand (time(NULL));
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
