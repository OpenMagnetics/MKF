#include "NumberTurns.h"
#include "TestingUtils.h"

#include <UnitTest++.h>
#include <vector>

SUITE(NumberTurns) {
    TEST(Number_Turns_Inductor) {
        OpenMagnetics::DesignRequirements designRequirements;
        designRequirements.set_turns_ratios(std::vector<OpenMagnetics::DimensionWithTolerance>{});
        uint64_t initialPrimaryNumberTurns = 42;

        OpenMagnetics::NumberTurns numberTurns(initialPrimaryNumberTurns, designRequirements);
        std::vector<uint64_t> numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        CHECK(numberTurnsCombination[0] == initialPrimaryNumberTurns);
        numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        CHECK(numberTurnsCombination[0] == initialPrimaryNumberTurns + 1);
    }

    TEST(Number_Turns_Two_Windings_Turns_Ratio_1) {
        OpenMagnetics::DesignRequirements designRequirements;
        OpenMagnetics::DimensionWithTolerance turnsRatio;
        double turnsRatioValue = 1;
        turnsRatio.set_nominal(turnsRatioValue);
        designRequirements.set_turns_ratios(std::vector<OpenMagnetics::DimensionWithTolerance>{turnsRatio});
        uint64_t initialPrimaryNumberTurns = 42;

        OpenMagnetics::NumberTurns numberTurns(initialPrimaryNumberTurns, designRequirements);
        std::vector<uint64_t> numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        CHECK(numberTurnsCombination[0] == initialPrimaryNumberTurns);
        CHECK(numberTurnsCombination[1] == initialPrimaryNumberTurns * 1);
        numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        CHECK(numberTurnsCombination[0] == initialPrimaryNumberTurns + 1);
        CHECK(numberTurnsCombination[1] == (initialPrimaryNumberTurns + 1) * 1);
    }

    TEST(Number_Turns_Two_Windings_Turns_Ratio_8) {
        OpenMagnetics::DesignRequirements designRequirements;
        OpenMagnetics::DimensionWithTolerance turnsRatio;
        double turnsRatioValue = 8;
        turnsRatio.set_nominal(turnsRatioValue);
        turnsRatio.set_minimum(turnsRatioValue * 0.8);
        turnsRatio.set_maximum(turnsRatioValue * 1.2);
        designRequirements.set_turns_ratios(std::vector<OpenMagnetics::DimensionWithTolerance>{turnsRatio});
        uint64_t initialPrimaryNumberTurns = 42;

        OpenMagnetics::NumberTurns numberTurns(initialPrimaryNumberTurns, designRequirements);
        std::vector<uint64_t> numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        CHECK(numberTurnsCombination[0] == initialPrimaryNumberTurns);
        CHECK(OpenMagnetics::check_requirement(turnsRatio, double(numberTurnsCombination[0]) / numberTurnsCombination[1]));

        numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        CHECK(numberTurnsCombination[0] == initialPrimaryNumberTurns + 1);
        CHECK(OpenMagnetics::check_requirement(turnsRatio, double(numberTurnsCombination[0]) / numberTurnsCombination[1]));
    }

    TEST(Number_Turns_Two_Windings_Turns_Ratio_0_001) {
        OpenMagnetics::DesignRequirements designRequirements;
        OpenMagnetics::DimensionWithTolerance turnsRatio;
        double turnsRatioValue = 0.001;
        turnsRatio.set_nominal(turnsRatioValue);
        turnsRatio.set_minimum(turnsRatioValue * 0.8);
        turnsRatio.set_maximum(turnsRatioValue * 1.2);
        designRequirements.set_turns_ratios(std::vector<OpenMagnetics::DimensionWithTolerance>{turnsRatio});
        uint64_t initialPrimaryNumberTurns = 42;

        OpenMagnetics::NumberTurns numberTurns(initialPrimaryNumberTurns, designRequirements);
        std::vector<uint64_t> numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        CHECK(numberTurnsCombination[0] == initialPrimaryNumberTurns);
        CHECK(OpenMagnetics::check_requirement(turnsRatio, double(numberTurnsCombination[0]) / numberTurnsCombination[1]));

        numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        CHECK(numberTurnsCombination[0] == initialPrimaryNumberTurns + 1);
        CHECK(OpenMagnetics::check_requirement(turnsRatio, double(numberTurnsCombination[0]) / numberTurnsCombination[1]));
    }

    TEST(Number_Turns_Two_Windings_Turns_Ratio_Random) {
        srand (time(NULL));
        for (size_t i = 0; i < 1000; ++i)
        {
            OpenMagnetics::DesignRequirements designRequirements;
            OpenMagnetics::DimensionWithTolerance turnsRatio;
            double turnsRatioValue =  ((double) std::rand() / RAND_MAX) * (100 - 0.0001) + 0.0001;
            if (std::rand() % 2 == 0) {
                turnsRatioValue = 1 / turnsRatioValue;
            }
            turnsRatio.set_nominal(turnsRatioValue);
            turnsRatio.set_minimum(turnsRatioValue * 0.95);
            turnsRatio.set_maximum(turnsRatioValue * 1.05);
            designRequirements.set_turns_ratios(std::vector<OpenMagnetics::DimensionWithTolerance>{turnsRatio});
            uint64_t initialPrimaryNumberTurns = std::rand() % 100 + 1UL;
            
            OpenMagnetics::NumberTurns numberTurns(initialPrimaryNumberTurns, designRequirements);
            std::vector<uint64_t> numberTurnsCombination = numberTurns.get_next_number_turns_combination();
            CHECK(numberTurnsCombination[0] >= initialPrimaryNumberTurns);
            CHECK(OpenMagnetics::check_requirement(turnsRatio, double(numberTurnsCombination[0]) / numberTurnsCombination[1]));

            numberTurnsCombination = numberTurns.get_next_number_turns_combination();
            CHECK(OpenMagnetics::check_requirement(turnsRatio, double(numberTurnsCombination[0]) / numberTurnsCombination[1]));
        }
    }

    TEST(Number_Turns_Many_Windings_Turns_Ratio_Random) {
        srand (time(NULL));
        for (size_t i = 0; i < 1000; ++i)
        {
            OpenMagnetics::DesignRequirements designRequirements;
            std::vector<OpenMagnetics::DimensionWithTolerance> turnsRatios;
            size_t numberSecondaryWindings = std::rand() % 10;
            for (size_t turnRatioIndex = 0; turnRatioIndex < numberSecondaryWindings; ++turnRatioIndex) {
                OpenMagnetics::DimensionWithTolerance turnsRatio;
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
            
            OpenMagnetics::NumberTurns numberTurns(initialPrimaryNumberTurns, designRequirements);
            std::vector<uint64_t> numberTurnsCombination = numberTurns.get_next_number_turns_combination();
            CHECK(numberTurnsCombination[0] >= initialPrimaryNumberTurns);
            for (size_t turnRatioIndex = 0; turnRatioIndex < turnsRatios.size(); ++turnRatioIndex)
            {
                CHECK(OpenMagnetics::check_requirement(turnsRatios[turnRatioIndex], double(numberTurnsCombination[0]) / numberTurnsCombination[turnRatioIndex + 1]));
            }

            numberTurnsCombination = numberTurns.get_next_number_turns_combination();
            for (size_t turnRatioIndex = 0; turnRatioIndex < turnsRatios.size(); ++turnRatioIndex)
            {
                CHECK(OpenMagnetics::check_requirement(turnsRatios[turnRatioIndex], double(numberTurnsCombination[0]) / numberTurnsCombination[turnRatioIndex + 1]));
            }

        }
    }

    TEST(Number_Turns_Two_Windings_Turns_Ratio_Random_0) {
        srand (time(NULL));
        OpenMagnetics::DesignRequirements designRequirements;
        OpenMagnetics::DimensionWithTolerance turnsRatio;
        double turnsRatioValue = 78;
        turnsRatio.set_nominal(turnsRatioValue);
        turnsRatio.set_minimum(turnsRatioValue * 0.8);
        turnsRatio.set_maximum(turnsRatioValue * 1.2);
        designRequirements.set_turns_ratios(std::vector<OpenMagnetics::DimensionWithTolerance>{turnsRatio});
        uint64_t initialPrimaryNumberTurns = 40;
        
        OpenMagnetics::NumberTurns numberTurns(initialPrimaryNumberTurns, designRequirements);
        std::vector<uint64_t> numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        CHECK(numberTurnsCombination[0] >= initialPrimaryNumberTurns);
        CHECK(OpenMagnetics::check_requirement(turnsRatio, double(numberTurnsCombination[0]) / numberTurnsCombination[1]));

        numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        CHECK(OpenMagnetics::check_requirement(turnsRatio, double(numberTurnsCombination[0]) / numberTurnsCombination[1]));
    }

    TEST(Number_Turns_Two_Windings_Turns_Ratio_Random_1) {
        srand (time(NULL));
        OpenMagnetics::DesignRequirements designRequirements;
        OpenMagnetics::DimensionWithTolerance turnsRatio;
        double turnsRatioValue = 0.010101;
        turnsRatio.set_nominal(turnsRatioValue);
        turnsRatio.set_minimum(turnsRatioValue * 0.8);
        turnsRatio.set_maximum(turnsRatioValue * 1.2);
        designRequirements.set_turns_ratios(std::vector<OpenMagnetics::DimensionWithTolerance>{turnsRatio});
        uint64_t initialPrimaryNumberTurns = 60;
        
        OpenMagnetics::NumberTurns numberTurns(initialPrimaryNumberTurns, designRequirements);
        std::vector<uint64_t> numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        CHECK(numberTurnsCombination[0] >= initialPrimaryNumberTurns);
        CHECK(OpenMagnetics::check_requirement(turnsRatio, double(numberTurnsCombination[0]) / numberTurnsCombination[1]));

        numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        CHECK(OpenMagnetics::check_requirement(turnsRatio, double(numberTurnsCombination[0]) / numberTurnsCombination[1]));
    }

}