#include "constructive_models/NumberTurns.h"
#include "support/Settings.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <streambuf>
#include <vector>

namespace OpenMagnetics {

std::vector<uint64_t> NumberTurns::get_next_number_turns_combination(size_t multiple) {
    auto currentNumberTurns = _currentNumberTurns;
    increment_number_turns(multiple);

    return currentNumberTurns;
}

void NumberTurns::increment_number_turns(size_t multiple) {
    uint64_t timeout = 1000;
    uint64_t primaryNumberTurns = _currentNumberTurns[0];
    bool allRequirementsPassed = false;
    while (!allRequirementsPassed) {
        _currentNumberTurns.clear();
        allRequirementsPassed = true;
        primaryNumberTurns += multiple;
        while (primaryNumberTurns % multiple != 0) {
            primaryNumberTurns++;
        }

        _currentNumberTurns.push_back(primaryNumberTurns);
        for (size_t turnsRatioIndex = 0; turnsRatioIndex < _turnsRatios.size(); ++turnsRatioIndex) {
            auto turnsRatioRequirement = _turnsRatiosRequirements[turnsRatioIndex];
            auto turnsRatio = _turnsRatios[turnsRatioIndex];

            uint64_t numberTurns = round(primaryNumberTurns / turnsRatio);
            if (check_requirement(turnsRatioRequirement, double(primaryNumberTurns) / numberTurns)) {
                _currentNumberTurns.push_back(numberTurns);
            }
            else {
                allRequirementsPassed = false;
                break;
            }
        }
        timeout--;
        if (timeout == 0) {
            throw std::runtime_error("NumberTurns did not converge");
        }
    }
}

} // namespace OpenMagnetics
