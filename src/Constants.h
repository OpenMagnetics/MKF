#if !defined(CONSTANTS_H)
#    define CONSTANTS_H 1
#    pragma once

#    include <numbers>

namespace OpenMagnetics {
struct Constants {
    Constants() {};
    const double residualGap = 5e-6;
    const double minimumNonResidualGap = 0.1e-3;
    const double vacuum_permeability = 1.25663706212e-6;

    const double magnetic_flux_density_saturation = 0.352; // HARDCODED TODO: replace when materials are implemented

    const double spacer_protuding_percentage = 0.2;
    const double windingPainterScale = 1e6;

    const double number_points_samples_waveforms = 128;
    const double minimum_distributed_fringing_factor = 1.05;
    const double maximum_distributed_fringing_factor = 1.3;
    const double initial_gap_length_for_searching = 0.001;

    const double roshen_magnetic_field_strength_step = 0.1;
};
} // namespace OpenMagnetics

#endif