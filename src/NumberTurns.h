#pragma once
#include "MAS.hpp"
#include "Utils.h"
#include <streambuf>
#include <vector>

namespace OpenMagnetics {

class NumberTurns {
    private:
        std::vector<uint64_t> _currentNumberTurns;
        std::vector<DimensionWithTolerance> _turnsRatiosRequirements;
        std::vector<double> _turnsRatios;
    protected:
    public:
        NumberTurns(double initialPrimaryNumberTurns, DesignRequirements designRequirements) {
            _turnsRatiosRequirements = designRequirements.get_turns_ratios();
            for (auto& turnsRatioRequirement : _turnsRatiosRequirements) {
                auto turnsRatio = resolve_dimensional_values(turnsRatioRequirement, DimensionalValues::NOMINAL);
                if (turnsRatio <= 0) {
                    throw std::invalid_argument("turns ratio  must be greater than 0: " + std::to_string(turnsRatio));
                }
                _turnsRatios.push_back(turnsRatio);
            }
            uint64_t InitialPrimaryNumberTurnsMinusOne = initialPrimaryNumberTurns - 1; // Because looking afor a new one will increment turns by one
            _currentNumberTurns.push_back(InitialPrimaryNumberTurnsMinusOne);
            increment_number_turns();
        }

        NumberTurns(double initialPrimaryNumberTurns) {
            uint64_t InitialPrimaryNumberTurnsMinusOne = initialPrimaryNumberTurns - 1; // Because looking afor a new one will increment turns by one
            _currentNumberTurns.push_back(InitialPrimaryNumberTurnsMinusOne);
            increment_number_turns();
        }

        std::vector<uint64_t> get_next_number_turns_combination(size_t multiple = 1);
        void increment_number_turns(size_t multiple = 1);
};


} // namespace OpenMagnetics