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

    // ---- Converter defensible defaults (used when the optional input is unset) ----
    // These are documented design defaults / standard references, centralized
    // here instead of as scattered magic literals. They are NOT silent
    // fallbacks for missing *required* engineering data (efficiency, fsw, etc.,
    // which throw via Topology::require_input); they are conventional defaults
    // for genuinely-optional refinements, each matching the topology's JSON
    // schema default or a cited recommendation.
    const double resonantQualityFactorDefaultLlc = 0.4;   // LLC / CLLLC schema default
    const double resonantQualityFactorDefaultCllc = 0.3;  // CLLC: Infineon AN (Q in 0.2-0.4); matches cllcResonant.json
    const double resonantQualityFactorDefaultSrc = 2.0;   // SRC schema default
    const double commonModeChokeLineImpedanceDefault = 50.0;  // CISPR-16 / LISN reference impedance per line [Ohm]
    const double currentTransformerMinimumMagnetizingInductance = 1e-6;  // CT magnetizing-inductance design floor [H]
    const double coupledInductorCouplingCoefficientDefault = 0.999;  // near-unity coupling (Zeta etc.) when unspecified

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
    // Core coating thickness per surface, in m, used when a core's coating is given
    // only by name (epoxy/parylene) with no explicit thickness — the bulk of the
    // catalogue. Datasheet-sourced (see memory/core-coating-thicknesses.md):
    //   - parylene: "0.0005\" Parylene C" finish, conformal/size-independent
    //               (Micrometals datasheets, e.g. TR11-17/94; = Magnetics AY finish)
    //   - epoxy:    derived from Magnetics 2017 Powder Core Catalog p.14 — epoxy adds
    //               0.18 mm (0.007") more OD than parylene (0.089 mm/side) → ~0.10 mm
    const double defaultParyleneCoreCoatingThickness = 12.7e-6;  // 0.0005"
    const double defaultEpoxyCoreCoatingThickness = 0.1e-3;      // 0.10 mm, POWDER cores
    // Epoxy on a FERRITE ring core is about 3x the powder-core film, so it is its own value.
    // Both derivations below are tolerance-free -- they compare a coated limit against the
    // UNCOATED limit in the same direction, so the bare core's dimensional tolerance cancels
    // instead of being counted as coating (comparing a coated max against a bare NOMINAL does
    // not, and on a part toleranced at a few percent the tolerance swamps the film):
    //   - TDK R 50.0x30.0x20.0 (B64290L0082) gives uncoated 50.0+/-1.0, 30.0+/-0.7, 20.0+/-0.5
    //     against coated 51.8 max / 28.5 min / 21.3 max. Against the uncoated LIMITS that is
    //     (51.8-51.0)/2, (29.3-28.5)/2 and (21.3-20.5)/2 = 0.400 mm on all three axes, which is
    //     exactly the "< 0.4 mm" epoxy layer TDK states in its ring-core general information.
    //   - Ferroxcube DIMENSIONS the jacket on its TN ring-core drawings as "(~0,3)" mm.
    // 0.4 mm is the limit, 0.3 mm the drawn nominal, so take the nominal.
    const double defaultFerriteEpoxyCoreCoatingThickness = 0.3e-3;  // Ferroxcube TN drawings
    // Insulation material (in insulation_materials.ndjson) that supplies the relative
    // permittivity / dielectric strength for each core coating type given only by name.
    const std::string defaultEpoxyCoreCoatingMaterial = "epoxy";        // eps_r 3.6
    const std::string defaultParyleneCoreCoatingMaterial = "parylene";  // eps_r 3.1 (Parylene C)
    const std::string defaultNylonCoreCoatingMaterial = "Nylon 6.6";    // eps_r 3.4
    // When a toroid carries no coating data, the default coating TYPE is chosen by size,
    // mirroring manufacturer practice (Micrometals: toroids with OD <= 0.20" / 5.08 mm
    // are vacuum-coated with Parylene C, larger cores get a two-color epoxy finish — the
    // thin parylene film does not consume the small bore). Thickness stays a flat
    // per-type value; the buildup is not size-scaled (Fair-Rite/Magnetics catalogs).
    const double defaultToroidParyleneMaximumOuterDiameter = 5.08e-3;  // 0.20", POWDER cores
    // Ferrite ring cores cross over to epoxy at a different, larger size than powder cores do,
    // and getting this wrong picks a coating ~24x the right thickness for a 6 mm ferrite toroid.
    // TDK states the split explicitly (parylene for < R 9.53, epoxy for R 9.53 and above), and
    // Magnetics' 2022 ferrite catalogue confirms it: its toroid tables carry coating code Y
    // (parylene) up to 7.62 mm OD and code Z (epoxy) from 9.53 mm OD upwards.
    const double defaultFerriteToroidParyleneMaximumOuterDiameter = 9.53e-3;  // TDK "R 9.53"
    // Radius of a coated ring core's cross-section edge, in m: what a turn is pulled over on a
    // toroid. NOT the same datum as the coating thicknesses above, which are the dielectric path
    // normal to the flat faces; this is the curvature the jacket forms AT the corner, on top of a
    // ferrite edge that was already broken. Nobody dimensions the bare ferrite edge -- IEC
    // 62317-12 gives a ring core as A/B/C and admits the chamfer only through the effective
    // height, the MMPA/IMA toroid specification says only that the corners "shall not be sharp or
    // rough", and Fair-Rite says its toroids are "supplied burnished to break sharp edges" -- but
    // the jacket over it IS dimensioned: Ferroxcube draws PA11 at ~0,3 mm on its TN ring-core
    // datasheets, and TDK bounds ring-core epoxy at < 0,4 mm. Take the drawn value; the tumbled
    // or chamfered ferrite underneath only adds to it, so this stays a lower bound.
    // Core::get_toroid_edge_radius() clamps it to what the ring section can geometrically carry.
    const double defaultToroidRingEdgeRadius = 0.3e-3;  // Ferroxcube TN drawings, TDK < 0.4 mm
    // Small ring cores are NOT chamfered. TDK's ring-core general information gives the edge
    // treatment by size tier -- small: "edges rounded by tumbling", medium: "chamfer on edges
    // and/or radius on the surface", medium/big: "chamfer on edges" -- so the drawn chamfer above
    // does not describe a small core, and the parylene film such a core carries (12,7 um) is only
    // the jacket, not the broken ferrite edge under it.
    //
    // NOTHING publishes the tumbled radius. IEC 62317-12 dimensions no edge at all, the MMPA/IMA
    // toroid specification asks only that corners "shall not be sharp or rough", Fair-Rite says
    // its toroids are "supplied burnished to break sharp edges", and neither TDK's nor
    // Ferroxcube's part datasheets carry a callout. The only published figure found is from the
    // ferrite patent literature, which puts the corner curvature of a ferrite core in the range
    // 0,02 to 0,2 mm. This takes the LOW end of that range, which is the conservative choice: it
    // is the tightest edge a tumbled ferrite corner is reported to have.
    //
    // CAVEAT, because this one is weak and it gates hard: the range is quoted for ferrite core
    // LEGS rather than ring-core edges, and at 0,02 mm the windability gate admits only a 0,04 mm
    // conductor, which is finer than small toroids are actually wound with. Treat a rejection on
    // a tumbled core as "unverified", not as a physical verdict, until a vendor gives a number.
    const double defaultTumbledToroidEdgeRadius = 0.02e-3;  // patent literature, low end of range
    const double overlappingFactorSurroundingTurns = 0.7;
    const WireStandard commonWireStandard = WireStandard::NEMA_MW_1000_C;
    const WiringTechnology wiringTechnology = WiringTechnology::WOUND;
    const double pcbInsulationThickness = 100e-6;
    const double minimumWireToWireDistance = 90e-6;
    const double minimumBorderToWireDistance = 90e-6;
    const double coreToLayerDistance = 250e-6;
    // Default PCB fabrication class for planar coils the advisers wind without a MAS group.pcb (MAS-RFC 0012 makes
    // pcb mandatory on printed groups). Ordinary 2-layer-capable fab: 0.4/0.3 mm vias, 0.3 mm copper spacing.
    const double pcbViaDiameter = 0.4e-3;
    const double pcbViaDrillDiameter = 0.3e-3;
    const double pcbViaToVia = 0.3e-3;
    const double pcbViaToTrack = 0.3e-3;
    // Board outline of an adviser-proposed planar: the core footprint plus room for terminals on both sides.
    const double pcbTerminalAreaWidth = 10e-3;
    const double pcbEdgeMargin = 2e-3;

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

    // Hard ceiling on a gapped core's fringing factor for the POST-GAP
    // winding-killer guard (reject_winding_killing_gaps). The gap is finalized
    // after the early fringing pass (add_initial_turns_by_inductance can balloon
    // it to clear saturation), so a second, hard reject runs on the final gap.
    //
    // ABT #14: the original single 1.6 ceiling let tiny gapped cores with a
    // fringing factor in (1.3, 1.6) survive, and the magnetic-adviser's final
    // scoring then preferred those small/cheap cores even though their solid-wire
    // gap-fringing proximity loss was ~18-20 W (litz does not fit their tiny
    // winding window). Tightening the primary ceiling to 1.3 prunes those
    // catastrophic cores up front, so the surviving pool fills with sane P/PQ/RM
    // cores that take litz at ~3 W — without any scoring change.
    //
    // 1.3 can be too strict for legitimate high-energy designs that genuinely
    // need a large gap (measured fringing factors of legitimate small inductors
    // reach ~1.47). To avoid returning ZERO cores in that case, the guard relaxes
    // to the looser ceiling below ONLY when the strict pass would reject every
    // candidate (logged, never silent — see reject_winding_killing_gaps). The
    // looser 1.6 still rejects the genuinely catastrophic gaps (the field-reported
    // winding-killer sat at ~1.71).
    const double coreAdviserWindingKillingFringingFactorLimit = 1.3;
    const double coreAdviserWindingKillingFringingFactorLimitRelaxed = 1.6;

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
