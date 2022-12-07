#if !defined(CONSTANTS_H)
#define CONSTANTS_H 1
#pragma once

#include <numbers>

namespace OpenMagnetics {
	struct Constants {
	    Constants() {};
	    const double residualGap = 5e-6;
	    const double minimumNonResidualGap = 0.1e-3;
	    const double vacuum_permeability = std::numbers::pi * 4e-7;
	};
}

#endif