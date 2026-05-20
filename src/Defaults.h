#if !defined(DEFAULTS_H)
#    define DEFAULTS_H 1
#    pragma once

#include "Models.h"
#include "cstdint"
#include <MAS.hpp>

using namespace MAS;

namespace OpenMagnetics {
struct Defaults {
    Defaults() {};
    const OpenMagnetics::CoreLossesModels coreLossesModelDefault = OpenMagnetics::CoreLossesModels::IGSE;
    const OpenMagnetics::CoreTemperatureModels coreTemperatureModelDefault = OpenMagnetics::CoreTemperatureModels::MANIKTALA;
    const OpenMagnetics::ReluctanceModels reluctanceModelDefault = OpenMagnetics::ReluctanceModels::ZHANG;
    // ALBACH + ROSHEN selected based on comprehensive validation (Jan 2026):
    // Lowest overall error: 9.6% ± 8.4% across 87 tests covering ROUND, LITZ, RECTANGULAR, and FOIL wires
    // See MKF/docs/WindingLossesModelValidation.md for full methodology and results
    const OpenMagnetics::MagneticFieldStrengthModels magneticFieldStrengthModelDefault = OpenMagnetics::MagneticFieldStrengthModels::ALBACH;
    const OpenMagnetics::MagneticFieldStrengthFringingEffectModels magneticFieldStrengthFringingEffectModelDefault = OpenMagnetics::MagneticFieldStrengthFringingEffectModels::ROSHEN;
    const OpenMagnetics::CoreThermalResistanceModels coreThermalResistanceModelDefault = OpenMagnetics::CoreThermalResistanceModels::MANIKTALA;
    const double maximumProportionMagneticFluxDensitySaturation = 0.7;
    const double coreAdviserFrequencyReference = 100000;
    const double coreAdviserMagneticFluxDensityReference = 0.5;
    const double coreAdviserThresholdValidity = 0.9;
    const double coreAdviserMaximumCoreTemperature = 150;
    const double coreAdviserMaximumPercentagePowerCoreLosses = 0.05;
    const uint64_t coreAdviserMaximumMagneticsAfterFiltering = 400;  // Higher limit to capture more single-stack options
    const uint64_t coreAdviserMaximumNumberStacks = 4;
    const double maximumCurrentDensity = 7000000;
    const double maximumCurrentDensityPlanar = 2000000;
    const double maximumEffectiveCurrentDensity = 12000000;
    const int maximumNumberParallels = 5;
    const double magneticFluxDensitySaturation = 0.5;
    const double ferriteInitialPermeability = 2000;
    const double ferriteSaturationFluxDensity = 0.35;  // Typical ferrite Bsat ~350 mT (e.g., 3C95, N95)
    const double magnetizingInductanceThresholdValidity = 0.25;
    const double selfResonantFrequencyMargin = 0.25;  // Maximum operating frequency as fraction of SRF (e.g., 0.25 = f_op < 0.25 * SRF)
    const double harmonicAmplitudeThreshold = 0.05;
    const double ambientTemperature = 25;
    const double measurementFrequency = 10000;
    const double maximumFrequency = 100e6;
    const int magneticFieldMirroringDimension = 1;
    const double maximumCoilPattern = 6;
    const WindingOrientation defaultRoundWindowSectionsOrientation = WindingOrientation::CONTIGUOUS;
    const WindingOrientation defaultRectangularWindowSectionsOrientation = WindingOrientation::OVERLAPPING;
    const CoilAlignment defaultRoundWindowSectionsAlignment = CoilAlignment::SPREAD;
    const CoilAlignment defaultRectangularWindowSectionsAlignment = CoilAlignment::INNER_OR_TOP;
    const std::string defaultEnamelledInsulationMaterial = "Polyurethane 155";
    const std::string defaultInsulationMaterial = "ETFE";
    const std::string defaultLayerInsulationMaterial = "Kapton HN";
    const std::string defaultConductorMaterial = "copper";
    const std::string defaultPcbInsulationMaterial = "FR4";
    const std::string defaultBobbinMaterial = "PET";
    const double overlappingFactorSurroundingTurns = 0.7;
    const WireStandard commonWireStandard = WireStandard::NEMA_MW_1000_C;
    const WiringTechnology wiringTechnology = WiringTechnology::WOUND;
    const double pcbInsulationThickness = 100e-6;
    const double minimumWireToWireDistance = 90e-6;
    const double minimumBorderToWireDistance = 90e-6;
    const double coreToLayerDistance = 250e-6;

    // Phase 4 (Group A): adviser sweep/optimization iteration budgets and step
    // factors. Previously magic-numbered in-line at the call sites; collected
    // here so they can be re-tuned (or eventually exposed via settings) in one
    // place.
    //
    // MagneticFilter skin-effect-aware sweep: how many turn-count iterations
    // the "increase N, recompute losses" loop will explore before giving up.
    // Used in MagneticCoreFilterMinimumLosses and the equivalent in
    // MagneticCoreFilterMinimumImpedance.
    const size_t coreAdviserSkinEffectMaxIterations = 10;
    // Step size for the N sweep, expressed as a fraction of the starting
    // turn count. With ~10 iterations a step of 1 only covers N..N+10, which
    // is too narrow for larger designs; ~10% of N_start gives geometric-ish
    // coverage out to ~2× N_start.
    const double coreAdviserSkinEffectTurnsStepFactor = 0.1;
    // Maximum iterations for the CMC impedance-fit "jump N" loop in
    // MagneticCoreFilterMinimumImpedance. Distinct from the skin-effect
    // sweep above; kept generous (100) because it's the *outer* fit loop.
    const int coilAdviserCmcImpedanceMaxIterations = 100;

    // CoreAdviser material-evaluation fan-out: how many top candidates from
    // the sorted material list to evaluate per magnetic. Different defaults
    // for the two code paths preserved:
    //  - add_powder_materials / alternative-material rank uses 10 (broader
    //    exploration, since powder lookups are cheap).
    //  - add_ferrite_materials_by_losses / by_impedance use 2 (each candidate
    //    materialises into a full advisory run, so cost is quadratic).
    const size_t coreAdviserAlternativeMaterialsNumberToUse = 10;
    const size_t coreAdviserFanOutMaterialsNumberToUse = 2;
    // Reference B used when building the synthetic sinusoidal excitation for
    // material *ranking* (NOT for the design itself). Lower than the design
    // reference because we want to compare materials at a B where Steinmetz
    // and proprietary loss models tend to agree.
    const double coreAdviserMagneticFluxDensityReferenceAlternative = 0.18;
    // Geometric / physical limits during gap optimisation. The "max practical
    // gap as a fraction of column width" rejects designs whose gap would
    // require an unphysical core geometry. The "max fringing factor" is the
    // threshold beyond which gap-fringing flux is considered excessive and
    // penalised in the analytical cost function.
    const double coreAdviserMaxPracticalGapColumnWidthFraction = 0.5;
    const double coreAdviserMaxFringingFactor = 0.25;

    // Phase 4 (Group B): cross-referencer scoring-normalisation floors.
    //
    // crossReferencerScoringAbsoluteFloor:
    //   Hard lower bound on the per-filter "minimum scoring" used for
    //   log-/linear-normalising. Prevents log10(0) and division by ~0 when
    //   a filter genuinely scores some candidate at 0 or below the noise
    //   floor. Was the hardcoded literal 1e-10 throughout.
    //
    // crossReferencerScoringDataRelativeFloorRatio:
    //   Floor on the minimum scoring expressed as a fraction of the maximum
    //   scoring in the same filter. Acts as a dynamic-range cap: if the
    //   spread between min and max is wider than 1/ratio, the minimum is
    //   raised so normalisation doesn't compress everything into ~1.0.
    //   Was the hardcoded literal 1e-6.
    //
    // crossReferencerNeutralScoreWhenEqual:
    //   Neutral mid-range score awarded when min==max (every candidate
    //   scored identically on this filter, so the filter carries no
    //   discriminating information). Was the hardcoded literal 0.5.
    const double crossReferencerScoringAbsoluteFloor = 1e-10;
    const double crossReferencerScoringDataRelativeFloorRatio = 1e-6;
    const double crossReferencerNeutralScoreWhenEqual = 0.5;
};

// Phase 3 (F7): canonical name for the "shape-only pre-filter" sentinel used
// by CoreAdviser. A Core/Bobbin/Wire whose name is exactly this string is a
// placeholder that has been materialised at the shape level but does NOT yet
// have a real material/bobbin/wire bound. Real material is bound later at
// the fan-out (see CoreAdviser.cpp::add_ferrite_materials_by_*).
// Some filters need to treat the Dummy sentinel specially to avoid attempting
// material-dependent calculations on an unbound placeholder.
// Bare string literal "Dummy" must not be reintroduced - use this constant so
// a future rename can land in one place.
inline constexpr const char* DUMMY_SENTINEL_NAME = "Dummy";

} // namespace OpenMagnetics

#endif
