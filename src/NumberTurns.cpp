#include "NumberTurns.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>

namespace OpenMagnetics {

std::vector<uint64_t> NumberTurns::get_next_number_turns_combination() {
    auto currentNumberTurns = _currentNumberTurns;
    increment_number_turns();

    return currentNumberTurns;
}

void NumberTurns::increment_number_turns() {
    uint64_t primaryNumberTurns = _currentNumberTurns[0];
    bool allRequirementsPassed = false;
    while (!allRequirementsPassed) {
        _currentNumberTurns.clear();
        allRequirementsPassed = true;
        primaryNumberTurns++;
        _currentNumberTurns.push_back(primaryNumberTurns);
        for (auto& turnsRatioRequirement : _turnsRatios){
            auto turnsRatio = resolve_dimensional_values(turnsRatioRequirement, DimensionalValues::NOMINAL);
            uint64_t numberTurns = round(primaryNumberTurns / turnsRatio);
            // std::cout << "numberTurns: " << numberTurns << std::endl;
            // std::cout << "primaryNumberTurns: " << primaryNumberTurns << std::endl;
            // std::cout << "turnsRatio: " << turnsRatio << std::endl;

            if (check_requirement(turnsRatioRequirement, double(primaryNumberTurns) / numberTurns)) {
                _currentNumberTurns.push_back(numberTurns);
            }
            else {
                allRequirementsPassed = false;
                break;
            }
        }
    }
}

} // namespace OpenMagnetics
