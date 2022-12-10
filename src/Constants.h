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

		const double absolute_permeability = vacuum_permeability * 2000; // HARDCODED TODO: replace when materials are implemented
		const double magnetic_flux_density_saturation = 0.352; // HARDCODED TODO: replace when materials are implemented
	};
}

#endif