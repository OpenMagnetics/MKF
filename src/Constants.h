#if !defined(CONSTANTS_H)
#define CONSTANTS_H 1
#pragma once

#include <numbers>
#include <vector>
#include <string>


namespace OpenMagnetics {
struct Constants {
    Constants() {};
    const double residualGap = 5e-6;
    const double minimumNonResidualGap = 0.1e-3;
    const double vacuum_permeability = 1.25663706212e-6;

    const double magnetic_flux_density_saturation = 0.352; // HARDCODED TODO: replace when materials are implemented

    const double spacer_protuding_percentage = 0.2;
    const double windingPainterScale = 1e6;
    const std::vector<std::string> windingPainterColorsScaleSections = {"#539796",
                                                                        "#0e1919",
                                                                        "#71b1b0",
                                                                        "#1c3232",
                                                                        "#94c4c4",
                                                                        "#2a4c4b",
                                                                        "#b8d8d7",
                                                                        "#376564",
                                                                        "#dbebeb",
                                                                        "#457e7d"};
    const std::vector<std::string> windingPainterColorsScaleLayers = {"#B18AEA",
                                                                      "#1b0935",
                                                                      "#c1a1ee",
                                                                      "#361369",
                                                                      "#d0b9f2",
                                                                      "#511c9e",
                                                                      "#e0d0f7",
                                                                      "#6c26d2",
                                                                      "#efe8fb",
                                                                      "#8e55e1"};
    const std::vector<std::string> windingPainterColorsScaleTurns = {"#FFB94E",
                                                                     "#372200",
                                                                     "#ffc771",
                                                                     "#6f4300",
                                                                     "#ffd595",
                                                                     "#a76500",
                                                                     "#ffe3b8",
                                                                     "#de8600",
                                                                     "#fff1dc",
                                                                     "#ffa317"};

    const double number_points_samples_waveforms = 128;
    const double minimum_distributed_fringing_factor = 1.05;
    const double maximum_distributed_fringing_factor = 1.3;
    const double initial_gap_length_for_searching = 0.001;

    const double roshen_magnetic_field_strength_step = 0.1;
};
} // namespace OpenMagnetics

#endif