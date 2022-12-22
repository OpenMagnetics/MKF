#if !defined(CONSTANTS_H)
#define CONSTANTS_H 1
#pragma once

#include <numbers>

namespace OpenMagnetics {
	struct Constants {
	    Constants() {};
	    const double residualGap = 5e-6;
	    const double minimumNonResidualGap = 0.1e-3;
	    const double vacuum_permeability = 1.25663706212e-6;

		const double absolute_permeability = vacuum_permeability * 2000; // HARDCODED TODO: replace when materials are implemented
		const double magnetic_flux_density_saturation = 0.352; // HARDCODED TODO: replace when materials are implemented

		const double spacer_protuding_percentage = 0.2;
	};
}

#endif