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
    const double vacuumPermeability = 1.25663706212e-6;
    const double vacuumPermittivity =  8.8541878128e-12;

    const double quasiStaticFrequencyLimit = 100;

    const double spacerProtudingPercentage = 0.2;
    const double coilPainterScale = 1e6;
    const std::vector<std::string> coilPainterColorsScaleSections = {"#539796",
                                                                     "#0e1919",
                                                                     "#71b1b0",
                                                                     "#1c3232",
                                                                     "#94c4c4",
                                                                     "#2a4c4b",
                                                                     "#b8d8d7",
                                                                     "#376564",
                                                                     "#dbebeb",
                                                                     "#457e7d"};
    const std::vector<std::string> coilPainterColorsScaleLayers = {"#B18AEA",
                                                                   "#1b0935",
                                                                   "#c1a1ee",
                                                                   "#361369",
                                                                   "#d0b9f2",
                                                                   "#511c9e",
                                                                   "#e0d0f7",
                                                                   "#6c26d2",
                                                                   "#efe8fb",
                                                                   "#8e55e1"};
    const std::vector<std::string> coilPainterColorsScaleTurns = {"#FFB94E",
                                                                  "#372200",
                                                                  "#ffc771",
                                                                  "#6f4300",
                                                                  "#ffd595",
                                                                  "#a76500",
                                                                  "#ffe3b8",
                                                                  "#de8600",
                                                                  "#fff1dc",
                                                                  "#ffa317"};

    const size_t numberPointsSampledWaveforms = 128;
    const double minimumDistributedFringingFactor = 1.05;
    const double maximumDistributedFringingFactor = 1.3;
    const double initialGapLengthForSearching = 0.001;

    const double roshenMagneticFieldStrengthStep = 0.1;
    const double foilToSectionMargin = 0.05;
    const double planarToSectionMargin = 0.05;
};
} // namespace OpenMagnetics

#endif