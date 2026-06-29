#include "constructive_models/NumberTurns.h"
#include "support/Settings.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <streambuf>
#include <vector>
#include "support/Exceptions.h"

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
            // For a large step-down ratio the secondary is tiny (3-4 turns), so the
            // integer round() of a single fixed primary can overshoot the requested
            // ratio by 10-15% while STILL sitting inside check_requirement's loose
            // ±25% nominal-only validity band. A realized ratio that high forces a
            // resonant converter to boost far off-resonance (high circulating current
            // -> low efficiency, or beyond the tank's reach) — abt #62. When the
            // requirement is a bare nominal (no explicit min/max tolerance band, i.e.
            // the caller wants "this exact ratio"), additionally require the realized
            // ratio to be within a TIGHT relative tolerance: the search then advances
            // to a slightly higher primary whose larger secondary quantises the ratio
            // finely (more primary turns -> lower flux density, so never a saturation
            // penalty). An explicit min/max band is honoured verbatim (it encodes an
            // intentional, possibly wider, tolerance) — only the nominal-only case is
            // tightened beyond the generic validity band.
            if (numberTurns == 0) {
                // Primary still too small to realise this ratio at all; keep growing it
                // (and guard the division below against a zero secondary).
                allRequirementsPassed = false;
                break;
            }
            double realizedRatio = double(primaryNumberTurns) / numberTurns;
            bool ratioAcceptable = check_requirement(turnsRatioRequirement, realizedRatio);
            bool nominalOnly = !turnsRatioRequirement.get_minimum() &&
                               !turnsRatioRequirement.get_maximum() &&
                               turnsRatioRequirement.get_nominal();
            if (ratioAcceptable && nominalOnly) {
                ratioAcceptable = std::abs(realizedRatio - turnsRatio) <= _turnsRatioRelativeTolerance * turnsRatio;
            }
            if (ratioAcceptable) {
                _currentNumberTurns.push_back(numberTurns);
            }
            else {
                allRequirementsPassed = false;
                break;
            }
        }
        timeout--;
        if (timeout == 0) {
            throw CalculationException(ErrorCode::CALCULATION_DIVERGENCE, "NumberTurns did not converge");
        }
    }
}

} // namespace OpenMagnetics
