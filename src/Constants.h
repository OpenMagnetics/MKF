#if !defined(CONSTANTS_H)
#define CONSTANTS_H 1
#pragma once

#include <numbers>
#include <vector>
#include <string>

namespace OpenMagnetics {
struct Constants {
    Constants() {};
    ~Constants() {};
    const double residualGap = 5e-6;
    const double minimumNonResidualGap = 0.1e-3;
    const double vacuumPermeability = 1.25663706212e-6;
    const double vacuumPermittivity =  8.8541878128e-12;
    
    // Thermal constants
    const double gravityAcceleration = 9.80665;                    // m/s² - standard gravity
    const double kelvinOffset = 273.15;                            // K - offset from Celsius to Kelvin
    const double stefanBoltzmannConstant = 5.670374419e-8;         // W/(m²·K⁴) - Stefan-Boltzmann constant

    const double quasiStaticFrequencyLimit = 100;

    // Numerical tolerances
    const double defaultConvergenceTolerance = 0.001;
    const double strictConvergenceTolerance = 1e-6;
    const double looseConvergenceTolerance = 0.01;
    const double coordinateTolerance = 0.05;
    
    // Levenberg-Marquardt algorithm parameters
    const double lmInitMu = 1e-03;
    const double lmStopThreshold = 1e-25;
    const double lmDiffDelta = 1e-19;

    const double spacerProtudingPercentage = 0.2;
    const double coilPainterScale = 30000;
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

// ============================================================================
// IMP-2: Thermal analysis named constants (replaces magic numbers)
// ============================================================================
namespace ThermalDefaults {
    // Wire defaults
    constexpr double kWire_DefaultWidth = 0.001;           // [m] 1mm
    constexpr double kWire_DefaultHeight = 0.001;          // [m] 1mm
    constexpr double kWire_CopperThermalConductivity = 385.0;  // [W/(m·K)]
    constexpr double kWire_DefaultEnamelThickness = 0.00005;   // [m] 50µm
    constexpr double kWire_DefaultEnamelConductivity = 0.2;    // [W/(m·K)]
    
    // Insulation & Bobbin
    constexpr double kBobbin_ThermalConductivity = 0.2;    // [W/(m·K)] nylon/PBT
    constexpr double kInsulation_DefaultConductivity = 0.2;// [W/(m·K)] polyimide
    constexpr double kInsulation_DefaultThickness = 0.0001;// [m] 0.1mm
    
    // Core
    constexpr double kCore_FerriteThermalConductivity = 4.0;   // [W/(m·K)] 3C90/3C95
    
    // Convection & Radiation
    constexpr double kConvection_DefaultEmissivity = 0.9;      // [-] dark/matte
    constexpr double kConvection_InitialDeltaT = 30.0;         // [°C] initial estimate
    constexpr double kConvection_MinNaturalH = 2.0;            // [W/(m²·K)] IMP-7
    constexpr double kConvection_MinForcedH = 10.0;            // [W/(m²·K)]
    
    // TIM
    constexpr double kTIM_DefaultResistance = 0.5;             // [K/W]
    
    // Solver
    constexpr double kSolver_MinConductance = 1e-9;            // [W/K]
} // namespace ThermalDefaults


// Input processing constants
namespace InputConstants {
    constexpr double kPowerConductionThreshold = 0.01;
    constexpr double kVoltageCloseToZeroThreshold = 0.05;
    constexpr double kZeroCrossingPointRatio = 0.02;
    constexpr double kNonConductingPointRatio = 0.1;
    constexpr double kSlopeComparisonThreshold = 0.01;
    constexpr double kSinusoidalFitErrorThreshold = 0.05;
    constexpr double kDutyCycleLowerBound = 0.03;
    constexpr double kDutyCycleUpperBound = 0.97;
    constexpr double kDutyCycleOnThreshold = 0.05;
}

} // namespace OpenMagnetics

#endif