#pragma once
#include "MAS.hpp"
#include "Utils.h"
#include <streambuf>
#include <vector>

namespace OpenMagnetics {

class NumberTurns {
    private:
        std::vector<int64_t> _currentNumberTurns;
        std::vector<DimensionWithTolerance> _turnsRatios;
    protected:
    public:
        NumberTurns(double initialPrimaryNumberTurns, DesignRequirements designRequirements) {
            _turnsRatios = designRequirements.get_turns_ratios();
            int64_t InitialPrimaryNumberTurnsMinusOne = initialPrimaryNumberTurns - 1; // Because looking afor a new one will increment turns by one
            _currentNumberTurns.push_back(InitialPrimaryNumberTurnsMinusOne);
            increment_number_turns();
        }

        NumberTurns(double initialPrimaryNumberTurns) {
            int64_t InitialPrimaryNumberTurnsMinusOne = initialPrimaryNumberTurns - 1; // Because looking afor a new one will increment turns by one
            _currentNumberTurns.push_back(InitialPrimaryNumberTurnsMinusOne);
            increment_number_turns();
        }

        std::vector<int64_t> get_next_number_turns_combination();
        void increment_number_turns();
};


} // namespace OpenMagnetics