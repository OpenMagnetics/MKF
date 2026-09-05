#include "physical_models/MagnetizingInductance.h"

#include "processors/Inputs.h"
#include "physical_models/ReluctanceNetwork.h"
#include "physical_models/MagneticField.h"
#include "physical_models/Reluctance.h"
#include "physical_models/InitialPermeability.h"
#include "constructive_models/CorePiece.h"
#include "support/Settings.h"
#include "support/Utils.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <streambuf>
#include <vector>
#include "support/Exceptions.h"

namespace OpenMagnetics {

double calculate_air_inductance(int64_t numberTurnsPrimary, Core core) {

    auto bobbin = Bobbin::create_quick_bobbin(core);
    auto bobbinColumnDepth = bobbin.get_processed_description()->get_column_depth();
    auto bobbinColumnWidth = bobbin.get_processed_description()->get_column_width().value();
    auto bobbinWindingWindowWidth = bobbin.get_winding_window_dimensions()[0];
    auto bobbinWindingWindowHeight = bobbin.get_winding_window_dimensions()[1];
    auto meanLengthRadius = (bobbinColumnDepth + bobbinColumnWidth) / 2 + bobbinWindingWindowWidth / 4;

    double coilInternalArea = std::numbers::pi * pow(meanLengthRadius, 2);

    double coreColumnArea = core.get_processed_description()->get_columns()[0].get_area();

    double airAreaProportion = (coilInternalArea - coreColumnArea) / coilInternalArea;

    return Constants().vacuumPermeability * pow(numberTurnsPrimary, 2) * (coilInternalArea * airAreaProportion * 2) / bobbinWindingWindowHeight;
}

std::pair<MagnetizingInductanceOutput, SignalDescriptor> MagnetizingInductance::calculate_inductance_and_magnetic_flux_density(Magnetic magnetic, OperatingPoint* operatingPoint) {
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();
    return calculate_inductance_and_magnetic_flux_density(core, coil, operatingPoint);
}



// ===================== Open-core (drum / rod) magnetizing inductance, ABT #331 =====================
// An open shape's magnetic circuit closes through the surrounding air, so the closed-circuit
// reluctance path (le/Ae of a closed core) does not apply — no closed le exists. The literature
// model is the demagnetising-factor effective permeability (Bozorth, "Ferromagnetism", 1945;
// cylinder refinements in Chen/Brug/Goldfarb, IEEE Trans. Magn. 27 (1991) 3601):
//
//     mu_rod = mu_i / (1 + N_d (mu_i - 1)),   N_d = axial factor of the equivalent spheroid.
//
// For a DRUM (post + two flanges) no exact closed form exists, so the estimate is the log-midpoint
// of two PHYSICAL BOUNDS:
//   upper: the flange-envelope spheroid (treats the whole drum as solid ferrite of flange
//          diameter) — measured bias +14% on the validation set;
//   lower: series reluctance of the envelope air return plus the ferrite post — bias -12%.
// The midpoint validates at mean 5.4% / max 9.9% against the published AL of the four Fair-Rite
// bobbins whose dimension mapping is weight-verified (9643001015: 38 nH, 9677282509: 95 nH,
// 9677182209: 65 nH, 9677282009: 100 nH) — see the pinned test. A key property this reproduces:
// AL is nearly material-independent (Fair-Rite publishes the SAME AL for 43 (mu_i 800) and 77
// (mu_i 2000) variants of one geometry) because for mu_i >> 1/N_d the result saturates at the
// geometry-set limit ~1/N_d.
static double open_core_axial_demagnetizing_factor(double aspectLengthOverDiameter) {
    double m = aspectLengthOverDiameter;
    if (fabs(m - 1) < 1e-4) {
        return 1.0 / 3;
    }
    if (m > 1) {
        double s = sqrt(m * m - 1);
        return (1 / (m * m - 1)) * ((m / s) * log(m + s) - 1);
    }
    double s = sqrt(1 - m * m);
    return (1 / (1 - m * m)) * (1 - m * acos(m) / s);
}

double MagnetizingInductance::calculate_open_core_magnetizing_inductance(Core core, double numberTurns, double temperature) {
    auto dimensions = flatten_dimensions(core.resolve_shape().get_dimensions().value());
    for (auto required : {"A", "B", "C", "D", "E", "F"}) {
        if (dimensions.find(required) == dimensions.end()) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                std::string("Open-core (drum) shape is missing dimension ") + required);
        }
    }
    // Only REAL gaps are rejected: process_gap() distributes residual bookkeeping entries onto
    // the columns of every non-toroidal core, drums included, so emptiness is not the test.
    for (auto& gap : core.get_functional_description().get_gapping()) {
        if (gap.get_type() != GapType::RESIDUAL && gap.get_length() > 1e-6) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                "An open-circuit (drum) core cannot be gapped: its return path is already air");
        }
    }
    // Asymmetric drums (A2 = second flange OD): the envelope the air return sees is between the
    // two flange discs; the geometric mean keeps the symmetric case exact and matched the WE-TI
    // reconstruction (both bounds are envelope-driven, so this is the sensitive choice).
    double flangeDiameter = dimensions["A"];
    if (dimensions.find("A2") != dimensions.end() && dimensions["A2"] > 0) {
        flangeDiameter = sqrt(dimensions["A"] * dimensions["A2"]);
    }
    double height = dimensions["B"];
    double postDiameter = dimensions["C"];
    double bore = (dimensions.find("H") != dimensions.end()) ? dimensions["H"] : 0.0;
    double grooveHeight = dimensions["E"];
    double initialPermeability = InitialPermeability::get_initial_permeability(core.resolve_material(), temperature);

    double vacuumPermeability = Constants().vacuumPermeability;
    double demagnetizingFactor = open_core_axial_demagnetizing_factor(height / flangeDiameter);
    double envelopeArea = std::numbers::pi / 4 * pow(flangeDiameter, 2);
    double postArea = std::numbers::pi / 4 * (pow(postDiameter, 2) - pow(bore, 2));
    double rodPermeability = initialPermeability / (1 + demagnetizingFactor * (initialPermeability - 1));

    double upperBound = vacuumPermeability * rodPermeability * envelopeArea / grooveHeight * pow(numberTurns, 2);
    double airReluctance = demagnetizingFactor * height / (vacuumPermeability * envelopeArea);
    double ferriteReluctance = height / (vacuumPermeability * initialPermeability * postArea);
    double lowerBound = pow(numberTurns, 2) / (airReluctance + ferriteReluctance);

    return sqrt(upperBound * lowerBound);
}

// ===================== ROD magnetizing inductance (open circuit), ABT #933 =====================
// A bare cylinder has NO return limb: every line of flux that leaves one end travels back through
// the surrounding air. The closed-circuit machinery (le/Ae over one mu) has nothing to describe
// here — le does not exist — so a rod uses the same demagnetising-factor treatment the drum does
// (Bozorth, "Ferromagnetism", 1945; cylinder refinements in Chen/Brug/Goldfarb, IEEE Trans. Magn.
// 27 (1991) 3601):
//
//     mu_rod = mu_i / (1 + N_d (mu_i - 1)),   N_d = axial demagnetising factor of the ROD.
//
// N_d is a property of the ROD's own slenderness B/A and of nothing else. That is the physical
// content of an open circuit and it is worth stating plainly, because it bounds what this model
// can and cannot explain: lengthening the WINDING on a fixed rod does not change N_d.
//
// The winding length enters through the other term. The rod-core inductor is the air solenoid
// lifted by mu_rod, and an air solenoid's inductance goes as N^2 A / l_winding, so a SHORT coil
// on a long rod has more inductance per turn than a full-length one. That is the only route by
// which the coil's own geometry reaches the answer, and it is bounded: over the full rod the two
// bracket bounds below collapse onto each other and the coil length drops out entirely.
//
// Same log-midpoint bracket idiom as the drum:
//   upper: the solenoid formula on the rod cross-section at mu_rod — exact for a long thin rod
//          fully wound, optimistic for a short coil because it ignores the flux that leaves the
//          rod alongside the winding;
//   lower: series reluctance of the external air return (N_d B / mu0 A_rod) and the rod itself
//          (B / mu0 mu_i A_rod) — pessimistic for the same reason, mirrored.
// For a fully wound rod the two agree to a fraction of a percent (for N_d = 0.2, mu_i = 400 they
// differ by 0.25%), so the bracket is not doing hidden work: it only opens up as the winding is
// shortened, which is exactly where the uncertainty actually lives.
static std::pair<double, double> rod_open_core_bounds(Core core, double numberTurns,
                                                     double windingLength, double temperature) {
    auto dimensions = flatten_dimensions(core.resolve_shape().get_dimensions().value());
    for (auto required : {"A", "B"}) {
        if (dimensions.find(required) == dimensions.end()) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                std::string("Open-core (rod) shape is missing dimension ") + required +
                " (a rod is A = diameter, B = length, H = optional bore)");
        }
    }
    // Only REAL gaps are rejected: process_gap() distributes residual bookkeeping entries onto
    // the columns of every non-toroidal core, rods included, so emptiness is not the test.
    for (auto& gap : core.get_functional_description().get_gapping()) {
        if (gap.get_type() != GapType::RESIDUAL && gap.get_length() > 1e-6) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                "An open-circuit (rod) core cannot be gapped: its return path is already air");
        }
    }
    double rodDiameter = dimensions["A"];
    double rodLength = dimensions["B"];
    double bore = (dimensions.find("H") != dimensions.end()) ? dimensions["H"] : 0.0;
    if (rodDiameter <= 0 || rodLength <= 0) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
            "Rod core has a non-positive diameter A or length B");
    }
    if (bore < 0 || bore >= rodDiameter) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
            "Rod core bore H must be non-negative and smaller than the diameter A");
    }
    if (!(windingLength > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
            "Rod magnetizing inductance needs the axial length of the winding, and none was "
            "available; a rod imposes no groove, so the coil must state it (bobbin winding "
            "window height)");
    }
    if (windingLength > rodLength * (1 + 1e-9)) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
            "The winding is longer (" + std::to_string(windingLength) + " m) than the rod it sits on (" +
            std::to_string(rodLength) + " m); the overhanging turns are air-cored and this model "
            "does not describe them");
    }
    double initialPermeability = InitialPermeability::get_initial_permeability(core.resolve_material(), temperature);

    double vacuumPermeability = Constants().vacuumPermeability;
    double demagnetizingFactor = open_core_axial_demagnetizing_factor(rodLength / rodDiameter);
    double rodArea = std::numbers::pi / 4 * (pow(rodDiameter, 2) - pow(bore, 2));
    double rodPermeability = initialPermeability / (1 + demagnetizingFactor * (initialPermeability - 1));

    double upperBound = vacuumPermeability * rodPermeability * rodArea / windingLength * pow(numberTurns, 2);
    double airReluctance = demagnetizingFactor * rodLength / (vacuumPermeability * rodArea);
    double ferriteReluctance = rodLength / (vacuumPermeability * initialPermeability * rodArea);
    double lowerBound = pow(numberTurns, 2) / (airReluctance + ferriteReluctance);

    return {lowerBound, upperBound};
}

double MagnetizingInductance::calculate_rod_magnetizing_inductance(Core core, double numberTurns,
                                                                  double windingLength, double temperature) {
    auto [lowerBound, upperBound] = rod_open_core_bounds(core, numberTurns, windingLength, temperature);
    return sqrt(upperBound * lowerBound);
}

// ABT #362: semi-shielded drum — a ferrite drum closed by a MAGNETIC-EPOXY shell. The circuit
// crosses TWO materials, so neither the closed-circuit path (one mu over le/Ae) nor a gap
// model fits: reluctance is applied PER SECTION using the piece's c1 split. The shell material
// rides the core coating {type: magneticEpoxy, material: <core material NAME>}; per the
// no-fallbacks rule every missing link throws — the model never assumes air or guesses a mu.
double MagnetizingInductance::calculate_semishielded_drum_magnetizing_inductance(Core core, double numberTurns, double temperature) {
    auto corePiece = CorePiece::factory(core.resolve_shape());
    auto mixedConstants = corePiece->get_mixed_material_constants();
    if (!mixedConstants) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
            "Semi-shielded drum piece did not expose its mixed-material shape constants");
    }
    double coreMaterialC1 = (*mixedConstants)[0];
    double shellMaterialC1 = (*mixedConstants)[2];
    double drumPermeability = InitialPermeability::get_initial_permeability(core.resolve_material(), temperature);

    auto coatingUnion = core.get_functional_description().get_coating();
    if (!coatingUnion) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
            "Semi-shielded drum requires a magneticEpoxy coating carrying the shell material — none present");
    }
    if (!std::holds_alternative<CoreCoating>(coatingUnion.value())) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
            "Semi-shielded drum coating must be an explicit CoreCoating object (type magneticEpoxy + "
            "shell material name), not a name-only coating");
    }
    auto coating = std::get<CoreCoating>(coatingUnion.value());
    if (!coating.get_type() || coating.get_type().value() != CoatingType::MAGNETIC_EPOXY) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
            "Semi-shielded drum requires coating type 'magneticEpoxy' (the glue shell IS the return path)");
    }
    if (!coating.get_material() || !std::holds_alternative<std::string>(coating.get_material().value())) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
            "Semi-shielded drum magneticEpoxy coating must reference the shell CORE material by name");
    }
    auto shellMaterial = find_core_material_by_name(std::get<std::string>(coating.get_material().value()));
    double shellPermeability = InitialPermeability::get_initial_permeability(shellMaterial, temperature);

    double vacuumPermeability = Constants().vacuumPermeability;
    double reluctance = (coreMaterialC1 / drumPermeability + shellMaterialC1 / shellPermeability) / vacuumPermeability;
    return pow(numberTurns, 2) / reluctance;
}

// ABT #576: a shielded drum and its closing ring are routinely DIFFERENT grades — WE-DPC-6040 is
// an ACME P47 MnZn drum (mu_i 3000) inside an ACME B45 NiZn ring (mu_i 450), a 6.7x difference
// across two material SYSTEMS. The grades are declared on functionalDescription.material as a
// LIST, primary piece first: ["P47", "B45"] is the drum then its ring. That field is the one the
// schema designates for analytical models; geometricalDescription is the CAD-facing view and is
// REGENERATED whenever absent, so a material declared only there would silently vanish and take
// the inductance 15% with it.
//
// A single material is not an error — the overwhelming majority of drumRing cores are one grade,
// and they keep the standard path. But a list that names MORE than two pieces is: a drumRing has
// exactly a drum and a ring, so anything else means the data does not describe this shape and we
// refuse rather than quietly using the first two.
static std::optional<CoreMaterial> resolve_drum_ring_ring_material(Core& core) {
    auto materials = core.resolve_materials();
    if (materials.size() == 1) {
        return std::nullopt;
    }
    if (materials.size() != 2) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
            "drumRing names " + std::to_string(materials.size()) + " materials; it has exactly two "
            "pieces (drum then ring), so the list must hold one or two entries");
    }
    return materials[1];
}

// Sectioned reluctance for a two-grade shielded drum: the drum sections take the drum's mu, the
// ring sections take the ring's mu, and the two STRUCTURAL annular clearances (ABT #366/#368,
// synthesised by Core::process_gap from A/K/D/F) are added on top exactly as before. Those
// clearances carry most of the reluctance in practice, which is why the single-material
// approximation is only ~15% out on inductance — but the grades differ far more than that for
// saturation and losses, so the split still matters.
double MagnetizingInductance::calculate_drum_ring_magnetizing_inductance(Core core,
                                                                        CoreMaterial ringMaterial,
                                                                        double numberTurns,
                                                                        double temperature) {
    auto corePiece = CorePiece::factory(core.resolve_shape());
    auto mixedConstants = corePiece->get_mixed_material_constants();
    if (!mixedConstants) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
            "drumRing piece did not expose its split drum/ring shape constants");
    }
    double drumC1 = (*mixedConstants)[0];
    double ringC1 = (*mixedConstants)[2];
    double drumPermeability = InitialPermeability::get_initial_permeability(core.resolve_material(), temperature);
    double ringPermeability = InitialPermeability::get_initial_permeability(ringMaterial, temperature);

    double vacuumPermeability = Constants().vacuumPermeability;
    double ferriteReluctance = (drumC1 / drumPermeability + ringC1 / ringPermeability) / vacuumPermeability;

    auto reluctanceModelForGaps = ReluctanceModel::factory(Defaults().reluctanceModelDefault);
    double clearanceReluctance = reluctanceModelForGaps->get_gapping_reluctance(core).get_gapping_reluctance().value();

    return pow(numberTurns, 2) / (ferriteReluctance + clearanceReluctance);
}

// ABT #1002: the WE moulded families are pressed in up to three steps and each pressing can be
// its own powder -- Inner / Outer in the MAPI and MAIA lists of parts, SUB / COR / COV (base,
// post, cover) in the MXGI one. The single fitted "effective" permeability the one-grade model
// absorbs those into is not a property of any powder, and it hides the fact that the post and
// the return path saturate at different currents. Here every region takes its own grade's mu at
// its own field: R = sum_r c1_r / (mu0 mu_r(H_r)), with H_r from the flux the DC magnetizing
// current drives through the whole series circuit, B_r = phi / Ae_r, iterated (damped, so a
// steep knee cannot make it ring) until every region's mu has settled.
double MagnetizingInductance::calculate_molded_magnetizing_inductance(Core core, double numberTurns, double temperature,
                                                                      std::optional<double> frequency,
                                                                      std::optional<double> magnetizingCurrentDcBias) {
    if (core.get_shape_family() != CoreShapeFamily::MOLDED) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
            "per-region moulded inductance asked of a core that is not a moulded body");
    }
    auto corePiece = CorePiece::factory(core.resolve_shape());
    auto regionConstants = corePiece->get_region_shape_constants();
    if (!regionConstants) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
            "moulded piece did not expose its per-region shape constants");
    }
    auto regionMaterials = core.resolve_region_materials();
    if (regionMaterials.size() == 1) {
        regionMaterials = std::vector<std::optional<CoreMaterial>>(regionConstants->size(), regionMaterials[0]);
    }
    if (regionMaterials.size() != regionConstants->size()) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
            "moulded body has " + std::to_string(regionConstants->size()) + " regions but its material "
            "list resolves to " + std::to_string(regionMaterials.size()));
    }
    if (!core.get_functional_description().get_gapping().empty()) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
            "a moulded body has no discrete gaps: its distributed gap lives in the powder");
    }

    double vacuumPermeability = Constants().vacuumPermeability;
    size_t numberRegions = regionConstants->size();
    auto permeabilityAt = [&](size_t regionIndex, std::optional<double> magneticFieldDcBias) {
        if (!regionMaterials[regionIndex]) {
            return 1.0;
        }
        return InitialPermeability::get_initial_permeability(regionMaterials[regionIndex].value(), temperature,
                                                             magneticFieldDcBias, frequency);
    };
    auto reluctanceOf = [&](const std::vector<double>& permeabilities) {
        double reluctance = 0;
        for (size_t regionIndex = 0; regionIndex < numberRegions; ++regionIndex) {
            reluctance += (*regionConstants)[regionIndex].c1 / (vacuumPermeability * permeabilities[regionIndex]);
        }
        return reluctance;
    };

    std::vector<double> permeabilities(numberRegions);
    for (size_t regionIndex = 0; regionIndex < numberRegions; ++regionIndex) {
        permeabilities[regionIndex] = permeabilityAt(regionIndex, std::nullopt);
    }
    double reluctance = reluctanceOf(permeabilities);

    if (magnetizingCurrentDcBias && *magnetizingCurrentDcBias != 0) {
        bool converged = false;
        for (size_t iteration = 0; iteration < 500 && !converged; ++iteration) {
            double magneticFlux = numberTurns * fabs(*magnetizingCurrentDcBias) / reluctance;
            double worstChange = 0;
            std::vector<double> updated(numberRegions);
            for (size_t regionIndex = 0; regionIndex < numberRegions; ++regionIndex) {
                auto& region = (*regionConstants)[regionIndex];
                double regionEffectiveArea = region.c1 / region.c2;
                double magneticFluxDensity = magneticFlux / regionEffectiveArea;
                double magneticFieldDcBias = magneticFluxDensity / (vacuumPermeability * permeabilities[regionIndex]);
                double target = permeabilityAt(regionIndex, magneticFieldDcBias);
                worstChange = std::max(worstChange, fabs(target - permeabilities[regionIndex]) / permeabilities[regionIndex]);
                // Geometric damping: mu(H) is decreasing, so the undamped fixed point can
                // oscillate around a steep knee; the half-step in log space cannot.
                updated[regionIndex] = sqrt(permeabilities[regionIndex] * target);
            }
            permeabilities = updated;
            reluctance = reluctanceOf(permeabilities);
            converged = worstChange < 1e-5;
        }
        if (!converged) {
            throw std::runtime_error("per-region moulded DC-bias iteration did not converge at " +
                                     std::to_string(*magnetizingCurrentDcBias) + " A");
        }
    }
    if (std::isnan(reluctance) || reluctance <= 0) {
        throw NaNResultException("moulded per-region reluctance must be a positive number");
    }
    return pow(numberTurns, 2) / reluctance;
}

// The DC component of the current that magnetizes the core, for the per-region moulded path:
// the magnetizing current when the caller already derived one, else the winding current of a
// single-winding part. Processed data is computed from the waveform when it is absent.
static std::optional<double> excitation_dc_current(const OperatingPointExcitation& excitation) {
    std::optional<SignalDescriptor> signal;
    if (excitation.get_magnetizing_current()) {
        signal = excitation.get_magnetizing_current();
    }
    else if (excitation.get_current()) {
        signal = excitation.get_current();
    }
    if (!signal) {
        return std::nullopt;
    }
    if (signal->get_processed()) {
        return signal->get_processed()->get_offset();
    }
    if (signal->get_waveform()) {
        // Same construction as the standard path: a power-of-two resample, its harmonics on the
        // descriptor, then the processed data -- calculate_processed_data reads the harmonics.
        auto sampled = Inputs::calculate_sampled_waveform(signal->get_waveform().value(), excitation.get_frequency());
        signal->set_harmonics(Inputs::calculate_harmonics_data(sampled, excitation.get_frequency()));
        auto processed = Inputs::calculate_processed_data(*signal, sampled, false);
        return processed.get_offset();
    }
    return std::nullopt;
}

// ABT #362/#331: the drum-family paths below return early with their own inductance model, so
// they must ALSO produce the magnetic flux density this function is contracted to return —
// MagneticSimulator feeds result.second straight into the core-loss stage, and a default-
// constructed SignalDescriptor makes IGSE throw bad_optional_access deep inside
// Inputs::get_magnetic_flux_density_peak_to_peak, naming nothing. Same construction as the
// main path: flux from the magnetizing current through the driving-point reluctance, then
// divided by the flux-carrying area.
static SignalDescriptor calculate_flux_density_for_family_model(double magnetizingInductance,
                                                                double numberTurns,
                                                                double fluxCarryingArea,
                                                                OperatingPoint* operatingPoint) {
    SignalDescriptor magneticFluxDensity;
    if (!operatingPoint || operatingPoint->get_mutable_excitations_per_winding().empty() ||
        !operatingPoint->get_mutable_excitations_per_winding()[0].get_magnetizing_current()) {
        return magneticFluxDensity;
    }
    double drivingPointReluctance = pow(numberTurns, 2) / magnetizingInductance;
    auto magneticFlux = OpenMagnetics::MagneticField::calculate_magnetic_flux(
        operatingPoint->get_mutable_excitations_per_winding()[0].get_magnetizing_current().value(),
        drivingPointReluctance, numberTurns);
    return OpenMagnetics::MagneticField::calculate_magnetic_flux_density(magneticFlux, fluxCarryingArea);
}

std::pair<MagnetizingInductanceOutput, SignalDescriptor> MagnetizingInductance::calculate_inductance_and_magnetic_flux_density(Core core, Coil coil, OperatingPoint* operatingPoint) {

    // ABT #417: a drumRing core's two structural annular-clearance gaps are DERIVED
    // (Core::process_gap synthesizes them from A/K/D/F — nothing is ever hand-authored
    // in functionalDescription.gapping for this family, unlike a normal gapped E-core).
    // Core's free from_json (used whenever a Core is deserialized as a MEMBER — e.g.
    // Magnetic::from_json, which every PyOM/WASM binding reaches) never calls
    // process_data()/process_gap(), unlike the Core(json) constructor. A drumRing core
    // arriving here via that path has gapping==[] and ReluctanceModel::get_gapping_reluctance
    // silently treats "empty" as "no gap" instead of "not yet derived", dropping the
    // dominant reluctance term and reporting an inductance 3.8-10.7x too high. gapping is
    // NEVER legitimately empty for this family (single- or dual-material), so self-heal
    // unconditionally here — process_gap() is idempotent and a no-op if already derived.
    if (core.get_shape_family() == CoreShapeFamily::DRUM_RING && core.get_functional_description().get_gapping().empty()) {
        core.process_gap_or_throw();
    }

    // Semi-shielded drums (ABT #362): mixed-material sectioned reluctance, drum mu + glue mu.
    if (core.get_shape_family() == CoreShapeFamily::DRUM_SEMISHIELDED) {
        double semishieldedTemperature = operatingPoint ? operatingPoint->get_conditions().get_ambient_temperature()
                                                        : Defaults().ambientTemperature;
        double numberTurnsSemishielded = coil.get_functional_description()[0].get_number_turns();
        double semishieldedInductance = calculate_semishielded_drum_magnetizing_inductance(core, numberTurnsSemishielded, semishieldedTemperature);
        MagnetizingInductanceOutput semishieldedOutput;
        DimensionWithTolerance semishieldedWithTolerance;
        semishieldedWithTolerance.set_nominal(semishieldedInductance);
        // Provisional envelope: unvalidated against the vendor set yet — the heimdall 290-part
        // sweep (ABT #362 acceptance) pins this; tighten when it lands.
        semishieldedWithTolerance.set_minimum(semishieldedInductance * 0.8);
        semishieldedWithTolerance.set_maximum(semishieldedInductance * 1.2);
        semishieldedOutput.set_magnetizing_inductance(semishieldedWithTolerance);
        semishieldedOutput.set_method_used("SemiShieldedMixedSectionReluctance");
        semishieldedOutput.set_origin(ResultOrigin::SIMULATION);
        std::pair<MagnetizingInductanceOutput, SignalDescriptor> semishieldedResult;
        semishieldedResult.first = semishieldedOutput;
        // Flux crosses the post: the drum's own central column carries it.
        semishieldedResult.second = calculate_flux_density_for_family_model(
            semishieldedInductance, numberTurnsSemishielded,
            core.get_columns()[0].get_area(), operatingPoint);
        return semishieldedResult;
    }

    // Shielded drums whose ring is a different grade from the drum (ABT #576). Gated on a
    // distinct ring grade actually being declared: with one grade this is a no-op and the core
    // falls through to the standard path below, so nothing already validated changes.
    if (core.get_shape_family() == CoreShapeFamily::DRUM_RING) {
        auto ringMaterial = resolve_drum_ring_ring_material(core);
        if (ringMaterial) {
            double drumRingTemperature = operatingPoint ? operatingPoint->get_conditions().get_ambient_temperature()
                                                        : Defaults().ambientTemperature;
            double numberTurnsDrumRing = coil.get_functional_description()[0].get_number_turns();
            double drumRingInductance = calculate_drum_ring_magnetizing_inductance(
                core, ringMaterial.value(), numberTurnsDrumRing, drumRingTemperature);
            MagnetizingInductanceOutput drumRingOutput;
            DimensionWithTolerance drumRingWithTolerance;
            drumRingWithTolerance.set_nominal(drumRingInductance);
            drumRingWithTolerance.set_minimum(drumRingInductance * 0.8);
            drumRingWithTolerance.set_maximum(drumRingInductance * 1.2);
            drumRingOutput.set_magnetizing_inductance(drumRingWithTolerance);
            drumRingOutput.set_method_used("DrumRingMixedSectionReluctance");
            drumRingOutput.set_origin(ResultOrigin::SIMULATION);
            std::pair<MagnetizingInductanceOutput, SignalDescriptor> drumRingResult;
            drumRingResult.first = drumRingOutput;
            // Flux crosses the post, same as every other drum-family model here.
            drumRingResult.second = calculate_flux_density_for_family_model(
                drumRingInductance, numberTurnsDrumRing,
                core.get_columns()[0].get_area(), operatingPoint);
            return drumRingResult;
        }
    }

    // Moulded bodies pressed from more than one powder (ABT #1002). Gated on the material list
    // actually resolving to more than one region, so every single-grade moulded core keeps the
    // standard path below unchanged.
    if (core.get_shape_family() == CoreShapeFamily::MOLDED && core.resolve_region_materials().size() > 1) {
        double moldedTemperature = operatingPoint ? operatingPoint->get_conditions().get_ambient_temperature()
                                                  : Defaults().ambientTemperature;
        double numberTurnsMolded = coil.get_functional_description()[0].get_number_turns();
        // The same reference frequency the standard path evaluates mu at when no operating point
        // is given, so a body that lists one grade three times reproduces the single-grade result.
        std::optional<double> moldedFrequency = Defaults().coreAdviserFrequencyReference;
        std::optional<double> moldedDcBias;
        if (operatingPoint && !operatingPoint->get_excitations_per_winding().empty()) {
            auto excitation = Inputs::get_primary_excitation(*operatingPoint);
            moldedFrequency = excitation.get_frequency();
            moldedDcBias = excitation_dc_current(excitation);
        }
        double moldedInductance = calculate_molded_magnetizing_inductance(
            core, numberTurnsMolded, moldedTemperature, moldedFrequency, moldedDcBias);
        MagnetizingInductanceOutput moldedOutput;
        DimensionWithTolerance moldedWithTolerance;
        moldedWithTolerance.set_nominal(moldedInductance);
        moldedWithTolerance.set_minimum(moldedInductance * 0.8);
        moldedWithTolerance.set_maximum(moldedInductance * 1.2);
        moldedOutput.set_magnetizing_inductance(moldedWithTolerance);
        moldedOutput.set_method_used("MoldedPerRegionReluctance");
        moldedOutput.set_origin(ResultOrigin::SIMULATION);
        std::pair<MagnetizingInductanceOutput, SignalDescriptor> moldedResult;
        moldedResult.first = moldedOutput;
        // Flux crosses the post, the region with the smallest section.
        moldedResult.second = calculate_flux_density_for_family_model(
            moldedInductance, numberTurnsMolded, core.get_columns()[0].get_area(), operatingPoint);
        return moldedResult;
    }

    // Rods (ABT #933): the most open shape there is — a bare cylinder with no return limb at all.
    // Its demagnetising model needs one thing a drum's does not, the axial length of the winding,
    // because a rod has no groove to imply it. That length is the coil's, so it is read here where
    // the coil is in scope and passed in; there is no default and no guess.
    if (core.get_shape_family() == CoreShapeFamily::ROD) {
        double rodTemperature = operatingPoint ? operatingPoint->get_conditions().get_ambient_temperature()
                                               : Defaults().ambientTemperature;
        double numberTurnsRod = coil.get_functional_description()[0].get_number_turns();
        double windingLength = coil.resolve_bobbin().get_winding_window_dimensions(0)[1];
        auto [rodLowerBound, rodUpperBound] = rod_open_core_bounds(core, numberTurnsRod, windingLength, rodTemperature);
        double rodInductance = sqrt(rodUpperBound * rodLowerBound);
        MagnetizingInductanceOutput rodOutput;
        DimensionWithTolerance rodWithTolerance;
        rodWithTolerance.set_nominal(rodInductance);
        // The two PHYSICAL bounds are themselves the honest uncertainty statement for this
        // geometry — no pinned percentage. Unlike the drum's validated +-12% this is not a
        // constant: it collapses to nothing on a fully wound rod and opens as the coil shortens,
        // which is exactly where the model's real uncertainty lives.
        rodWithTolerance.set_minimum(rodLowerBound);
        rodWithTolerance.set_maximum(rodUpperBound);
        rodOutput.set_magnetizing_inductance(rodWithTolerance);
        rodOutput.set_method_used("RodOpenCoreDemagnetizingFactor");
        rodOutput.set_origin(ResultOrigin::SIMULATION);
        std::pair<MagnetizingInductanceOutput, SignalDescriptor> rodResult;
        rodResult.first = rodOutput;
        // Flux crosses the rod itself (the bore is already subtracted from its column area).
        rodResult.second = calculate_flux_density_for_family_model(
            rodInductance, numberTurnsRod, core.get_columns()[0].get_area(), operatingPoint);
        return rodResult;
    }

    // Open shapes (drums, rods): route to the open-core model — the closed-circuit reluctance
    // machinery below would silently drop the dominant air-return reluctance (ABT #331).
    if (core.get_shape_family() == CoreShapeFamily::DRUM) {
        double openCoreTemperature = operatingPoint ? operatingPoint->get_conditions().get_ambient_temperature()
                                                    : Defaults().ambientTemperature;
        double numberTurnsOpenCore = coil.get_functional_description()[0].get_number_turns();
        double openCoreInductance = calculate_open_core_magnetizing_inductance(core, numberTurnsOpenCore, openCoreTemperature);
        MagnetizingInductanceOutput openCoreOutput;
        DimensionWithTolerance openCoreWithTolerance;
        openCoreWithTolerance.set_nominal(openCoreInductance);
        // Documented model envelope from the Fair-Rite validation set (max 9.9%).
        openCoreWithTolerance.set_minimum(openCoreInductance * 0.88);
        openCoreWithTolerance.set_maximum(openCoreInductance * 1.12);
        openCoreOutput.set_magnetizing_inductance(openCoreWithTolerance);
        openCoreOutput.set_method_used("OpenCoreDemagnetizingFactor");
        openCoreOutput.set_origin(ResultOrigin::SIMULATION);
        std::pair<MagnetizingInductanceOutput, SignalDescriptor> openCoreResult;
        openCoreResult.first = openCoreOutput;
        // Flux crosses the drum post (the bore is already subtracted from its column area).
        openCoreResult.second = calculate_flux_density_for_family_model(
            openCoreInductance, numberTurnsOpenCore,
            core.get_columns()[0].get_area(), operatingPoint);
        return openCoreResult;
    }


    double frequency = Defaults().coreAdviserFrequencyReference;
    double temperature = Defaults().ambientTemperature;

    if (operatingPoint) {

        temperature = operatingPoint->get_conditions().get_ambient_temperature();
        if (operatingPoint->get_mutable_excitations_per_winding().size() > 0) {
            frequency = operatingPoint->get_mutable_excitations_per_winding()[0].get_frequency();
            OperatingPointExcitation excitation = Inputs::get_primary_excitation(*operatingPoint);

            if (excitation.get_current()) {
                excitation.set_current(OpenMagnetics::Inputs::standardize_waveform(excitation.get_current().value(), frequency));
            }
            if (excitation.get_voltage()) {
                excitation.set_voltage(OpenMagnetics::Inputs::standardize_waveform(excitation.get_voltage().value(), frequency));
            }
            operatingPoint->get_mutable_excitations_per_winding()[0] = excitation;

            Inputs::make_waveform_size_power_of_two(operatingPoint);
        }
    }

    std::pair<MagnetizingInductanceOutput, SignalDescriptor> result;
    double numberWindings = coil.get_functional_description().size();
    double numberTurnsPrimary = coil.get_functional_description()[0].get_number_turns();
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    OpenMagnetics::InitialPermeability initialPermeability;
    double currentInitialPermeability;

    ReluctanceModels reluctanceModelEnum;
    from_json(_models["gapReluctance"], reluctanceModelEnum);
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(reluctanceModelEnum);
    double currentTotalReluctance;
    double modifiedTotalReluctance = 0;
    double modifiedMagnetizingInductance = 5e-3;
    double currentMagnetizingInductance;


    if (operatingPoint) {
        if (operatingPoint->get_mutable_excitations_per_winding().size() > 0) {
            OperatingPointExcitation excitation = Inputs::get_primary_excitation(*operatingPoint);
            if (!excitation.get_voltage()) {
                auto current = operatingPoint->get_mutable_excitations_per_winding()[0].get_current().value();
                auto currentWaveform = current.get_waveform().value();
                if (!is_size_power_of_2(currentWaveform.get_data())) {
                    auto currentSampledWaveform = Inputs::calculate_sampled_waveform(currentWaveform, frequency);
                    current.set_waveform(currentSampledWaveform);
                    operatingPoint->get_mutable_excitations_per_winding()[0].set_current(current);
                }
            }

            if (excitation.get_voltage()) {
                auto aux = operatingPoint->get_mutable_excitations_per_winding()[0].get_voltage().value().get_waveform().value();
                if (aux.get_data().size() > 0 && ((aux.get_data().size() & (aux.get_data().size() - 1)) != 0)) {
                    throw std::invalid_argument("voltage_data vector size is not a power of 2");
                }
            }
            if (excitation.get_current()) {
                auto aux = operatingPoint->get_mutable_excitations_per_winding()[0].get_current().value().get_waveform().value();
                if (aux.get_data().size() > 0 && ((aux.get_data().size() & (aux.get_data().size() - 1)) != 0)) {
                    throw std::invalid_argument("current_data vector size is not a power of 2");
                }
            }
            if (!excitation.get_voltage()) {
                if (!excitation.get_magnetizing_current()) {
                    Inputs::set_current_as_magnetizing_current(operatingPoint);
                }
                auto aux = operatingPoint->get_mutable_excitations_per_winding()[0].get_magnetizing_current().value().get_waveform().value();
                if (aux.get_data().size() > 0 && ((aux.get_data().size() & (aux.get_data().size() - 1)) != 0)) {
                    throw std::invalid_argument("magnetizing_current_data vector size from current is not a power of 2");
                }

            }
        }
    }

    MagnetizingInductanceOutput magnetizingInductanceOutput;
    currentInitialPermeability = initialPermeability.get_initial_permeability(core.resolve_material(), temperature, std::nullopt, frequency);
    magnetizingInductanceOutput = reluctanceModel->get_core_reluctance(core, currentInitialPermeability);
    modifiedTotalReluctance = magnetizingInductanceOutput.get_core_reluctance();
    modifiedMagnetizingInductance = pow(numberTurnsPrimary, 2) / modifiedTotalReluctance;

    size_t externalTimeout = 1;
    do {
        currentMagnetizingInductance = modifiedMagnetizingInductance;
        size_t internalTimeout = 1;

        do {
            currentTotalReluctance = modifiedTotalReluctance;
            modifiedMagnetizingInductance = pow(numberTurnsPrimary, 2) / currentTotalReluctance;


            if (operatingPoint) {
                if (operatingPoint->get_mutable_excitations_per_winding().size() > 0) {
                    OperatingPointExcitation excitation = Inputs::get_primary_excitation(*operatingPoint);

                    // If the converter model already computed the magnetizing current
                    // (e.g. DMC summing all winding currents), respect it — don't overwrite.
                    // BUT: only honor a preset MC when there is no voltage. When voltage is
                    // present, re-derive MC from voltage (V = L * dI/dt) using the runtime-
                    // computed magnetizing inductance — otherwise stale preset MCs (e.g.
                    // generated against a slightly different reluctance) cause B to drift
                    // and downstream loss tests to fail.
                    if (excitation.get_magnetizing_current() && !excitation.get_voltage()) {
                        // Already set — skip derivation. But upstream callers
                        // (MagneticField::get_magnetic_field_strength_gap and
                        // siblings) persist a *compressed* magnetizing_current
                        // (compress=true), which keeps only inflection points
                        // — typically a non-power-of-2 size like 23 for a
                        // triangular waveform. Resample to a dense
                        // power-of-2 waveform so the size-check gate below
                        // and the downstream FFT pipeline see a standardized
                        // contract.
                        auto presetMc = excitation.get_magnetizing_current().value();
                        if (presetMc.get_waveform()) {
                            auto presetWaveform = presetMc.get_waveform().value();
                            if (presetWaveform.get_data().size() > 0 && !is_size_power_of_2(presetWaveform.get_data())) {
                                if (!presetWaveform.get_time()) {
                                    auto stdMc = Inputs::standardize_waveform(presetMc, excitation.get_frequency());
                                    presetWaveform = stdMc.get_waveform().value();
                                }
                                auto sampled = Inputs::calculate_sampled_waveform(presetWaveform, excitation.get_frequency());
                                presetMc.set_waveform(sampled);
                                presetMc.set_harmonics(Inputs::calculate_harmonics_data(sampled, excitation.get_frequency()));
                                presetMc.set_processed(Inputs::calculate_processed_data(presetMc, sampled, false));
                                excitation.set_magnetizing_current(presetMc);
                                operatingPoint->get_mutable_excitations_per_winding()[0] = excitation;
                            }
                            else if (!presetMc.get_harmonics()) {
                                // Waveform is already power-of-2 but harmonics are missing (e.g. loaded from JSON
                                // with harmonics:null). Compute them so downstream frequency-domain paths work.
                                // Also recompute processed to avoid stale fields from the loaded JSON.
                                auto sampled = Inputs::calculate_sampled_waveform(presetWaveform, excitation.get_frequency());
                                presetMc.set_harmonics(Inputs::calculate_harmonics_data(sampled, excitation.get_frequency()));
                                presetMc.set_processed(Inputs::calculate_processed_data(presetMc, sampled, false));
                                excitation.set_magnetizing_current(presetMc);
                                operatingPoint->get_mutable_excitations_per_winding()[0] = excitation;
                            }
                        }
                    }
                    else if (numberWindings == 1 && excitation.get_current()) {
                        Inputs::set_current_as_magnetizing_current(operatingPoint);
                    }
                    // CMC check must come BEFORE is_multiport_inductor. CMCs
                    // have every winding on the same isolation side (L, N, PE
                    // all primary-side mains), which makes is_multiport_inductor
                    // return true — wrongly routing to the multiport path that
                    // uses the primary winding's raw current (DM + CM) as
                    // magnetizing current. That pumps the DM line current into
                    // the core flux calculation and overstates B by the ratio
                    // I_line / I_cm. The correct CMC path averages the winding
                    // currents and drops the DM DC offset.
                    else if (Inputs::can_be_common_mode_choke(*operatingPoint) && core.get_type() == CoreType::TOROIDAL) {
                        auto magnetizingCurrent = Inputs::get_common_mode_choke_magnetizing_current(*operatingPoint);
                        excitation.set_magnetizing_current(magnetizingCurrent);
                        operatingPoint->get_mutable_excitations_per_winding()[0] = excitation;
                    }
                    else if (Inputs::is_multiport_inductor(*operatingPoint, coil.get_isolation_sides())) {
                        auto magnetizingCurrent = Inputs::get_multiport_inductor_magnetizing_current(*operatingPoint);
                        excitation.set_magnetizing_current(magnetizingCurrent);
                        operatingPoint->get_mutable_excitations_per_winding()[0] = excitation;
                    }
                    else if (excitation.get_voltage()) {
                        auto voltage = operatingPoint->get_mutable_excitations_per_winding()[0].get_voltage().value();
                        auto sampledVoltageWaveform = Inputs::calculate_sampled_waveform(voltage.get_waveform().value(), frequency);

                        auto turnsRatios = coil.get_turns_ratios();
                        bool addOffset = Inputs::include_dc_offset_into_magnetizing_current(*operatingPoint, turnsRatios);

                        auto magnetizingCurrent = Inputs::calculate_magnetizing_current(excitation,
                                                                                                sampledVoltageWaveform,
                                                                                                modifiedMagnetizingInductance,
                                                                                                false,
                                                                                                addOffset,
                                                                                                operatingPoint->get_excitations_per_winding().size() > 1);

                        auto sampledMagnetizingCurrentWaveform = Inputs::calculate_sampled_waveform(magnetizingCurrent.get_waveform().value(), excitation.get_frequency());
                        // Replace the stored waveform with the resampled (power-of-2)
                        // version so the size-check gate below (and any downstream
                        // FFT pipeline) sees a standardized contract.
                        magnetizingCurrent.set_waveform(sampledMagnetizingCurrentWaveform);
                        magnetizingCurrent.set_harmonics(Inputs::calculate_harmonics_data(sampledMagnetizingCurrentWaveform, excitation.get_frequency()));
                        magnetizingCurrent.set_processed(Inputs::calculate_processed_data(magnetizingCurrent, sampledMagnetizingCurrentWaveform, false));

                        excitation.set_magnetizing_current(magnetizingCurrent);
                        operatingPoint->get_mutable_excitations_per_winding()[0] = excitation;
                    }

                    auto aux = operatingPoint->get_mutable_excitations_per_winding()[0].get_magnetizing_current().value().get_waveform().value();
                    if (aux.get_data().size() > 0 && ((aux.get_data().size() & (aux.get_data().size() - 1)) != 0)) {
                        throw std::invalid_argument("magnetizing_current_data vector size from voltage is not a power of 2 [size=" + std::to_string(aux.get_data().size()) + "]");
                    }

                    if (!operatingPoint->get_mutable_excitations_per_winding()[0].get_magnetizing_current()->get_waveform()->get_time()) {
                        auto magnetizingCurrent = Inputs::standardize_waveform(operatingPoint->get_mutable_excitations_per_winding()[0].get_magnetizing_current().value(), excitation.get_frequency());
                        operatingPoint->get_mutable_excitations_per_winding()[0].set_magnetizing_current(magnetizingCurrent);
                    }

                    auto magneticFlux = OpenMagnetics::MagneticField::calculate_magnetic_flux(operatingPoint->get_mutable_excitations_per_winding()[0].get_magnetizing_current().value(), currentTotalReluctance, numberTurnsPrimary);
                    auto magneticFluxDensity = OpenMagnetics::MagneticField::calculate_magnetic_flux_density(magneticFlux, effectiveArea);
                    result.second = magneticFluxDensity;
                    auto magneticFieldStrength = OpenMagnetics::MagneticField::calculate_magnetic_field_strength(magneticFluxDensity, currentInitialPermeability);
                    double switchingFrequency = Inputs::get_switching_frequency(operatingPoint->get_mutable_excitations_per_winding()[0]);

                    double hFieldDcBias = magneticFieldStrength.get_processed().value().get_offset();
                    if (!magneticFieldStrength.get_harmonics()) {
                        throw std::runtime_error("magneticFieldStrength has no harmonics — upstream magnetizing_current must provide a populated harmonics block (preset waveform missing harmonics?)");
                    }
                    if (magneticFieldStrength.get_harmonics().value().get_frequencies()[1] < switchingFrequency) {
                        for (size_t i = 0; i < magneticFieldStrength.get_harmonics().value().get_frequencies().size() - 1; ++i) {
                            if (magneticFieldStrength.get_harmonics().value().get_frequencies()[i] >= switchingFrequency) {
                                break;
                            }
                            hFieldDcBias = std::max(hFieldDcBias, magneticFieldStrength.get_harmonics().value().get_amplitudes()[i]);
                        }
                    }

                    currentInitialPermeability = initialPermeability.get_initial_permeability(core.resolve_material(), temperature, hFieldDcBias, frequency);

                    magnetizingInductanceOutput = reluctanceModel->get_core_reluctance(core, currentInitialPermeability);
                    modifiedTotalReluctance = magnetizingInductanceOutput.get_core_reluctance();
                    modifiedMagnetizingInductance = pow(numberTurnsPrimary, 2) / modifiedTotalReluctance;
                }
            }

            internalTimeout--;
            if (internalTimeout == 0) {
                break;
            }
        } while (fabs(currentTotalReluctance - modifiedTotalReluctance) / modifiedTotalReluctance >= 0.1);

        externalTimeout--;
        if (externalTimeout == 0) {
            break;
        }
    } while (fabs(currentMagnetizingInductance - modifiedMagnetizingInductance) / modifiedMagnetizingInductance >= 0.1);

    // Multi-column winding placement: the lumped N²/R model assumes the primary links
    // the main-column flux. When any winding is placed on another column, rebuild the
    // driving-point magnetizing inductance of the primary from the per-column
    // reluctance network at the converged permeability.
    {
        Magnetic magneticForPlacement;
        magneticForPlacement.set_core(core);
        magneticForPlacement.set_coil(coil);
        if (ReluctanceNetwork::has_non_main_placement(magneticForPlacement)) {
            ReluctanceNetwork magneticCircuit(core, magnetizingInductanceOutput.get_ungapped_core_reluctance().value(),
                                            magnetizingInductanceOutput.get_reluctance_per_gap().value_or(std::vector<AirGapReluctanceOutput>{}));
            auto inductanceMatrix = magneticCircuit.calculate_magnetizing_inductance_matrix(magneticForPlacement);
            modifiedMagnetizingInductance = inductanceMatrix[0][0];

            // Per-column flux density: the driven column carries the full
            // driving-point flux through ITS area. The lumped B = Φ/Ae both uses the
            // wrong reluctance for a lateral-driven leg and averages over the
            // effective area; saturation is checked in the actual column.
            if (operatingPoint && operatingPoint->get_mutable_excitations_per_winding().size() > 0 &&
                operatingPoint->get_mutable_excitations_per_winding()[0].get_magnetizing_current()) {
                auto columnIndexPerWinding = ReluctanceNetwork::resolve_winding_column_indexes(magneticForPlacement);
                double drivenColumnArea = core.get_columns()[columnIndexPerWinding[0]].get_area();
                double drivingPointReluctance = pow(numberTurnsPrimary, 2) / modifiedMagnetizingInductance;
                auto magneticFlux = OpenMagnetics::MagneticField::calculate_magnetic_flux(
                    operatingPoint->get_mutable_excitations_per_winding()[0].get_magnetizing_current().value(),
                    drivingPointReluctance, numberTurnsPrimary);
                result.second = OpenMagnetics::MagneticField::calculate_magnetic_flux_density(magneticFlux, drivenColumnArea);
            }
        }
    }

    if (operatingPoint) {
        if (operatingPoint->get_mutable_excitations_per_winding().size() > 0) {
            OperatingPointExcitation excitation = Inputs::get_primary_excitation(*operatingPoint);
            if (!excitation.get_voltage()) {
                operatingPoint->get_mutable_excitations_per_winding()[0].set_voltage(Inputs::calculate_induced_voltage(excitation, modifiedMagnetizingInductance));
            }
        }
    }

    auto& settings = Settings::GetInstance();
    if (settings.get_magnetizing_inductance_include_air_inductance()) {
        modifiedMagnetizingInductance += calculate_air_inductance(numberTurnsPrimary, core);
    }

    DimensionWithTolerance magnetizingInductanceWithTolerance;
    magnetizingInductanceWithTolerance.set_nominal(modifiedMagnetizingInductance);
    magnetizingInductanceOutput.set_magnetizing_inductance(magnetizingInductanceWithTolerance);

    result.first = magnetizingInductanceOutput;
    return result;
}

double MagnetizingInductance::calculate_inductance_air_solenoid(Magnetic magnetic) {
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();
    return calculate_inductance_air_solenoid(core, coil);
}


double MagnetizingInductance::calculate_inductance_air_solenoid(Core core, Coil coil) {
    double numberTurnsPrimary = coil.get_functional_description()[0].get_number_turns();

    ReluctanceModels reluctanceModelEnum;
    from_json(_models["gapReluctance"], reluctanceModelEnum);
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(reluctanceModelEnum);
    double airCoreReluctance = reluctanceModel->get_air_cored_reluctance(coil.resolve_bobbin());
    auto modifiedMagnetizingInductance = pow(numberTurnsPrimary, 2) / airCoreReluctance;
    return modifiedMagnetizingInductance;
}

MagnetizingInductanceOutput MagnetizingInductance::calculate_inductance_from_number_turns_and_gapping(Core core, Coil coil, OperatingPoint* operatingPoint) {
    auto inductance_and_magnetic_flux_density = calculate_inductance_and_magnetic_flux_density(core, coil, operatingPoint);

    return inductance_and_magnetic_flux_density.first;
}

MagnetizingInductanceOutput MagnetizingInductance::calculate_inductance_from_number_turns_and_gapping(Magnetic magnetic, OperatingPoint* operatingPoint) {
    auto inductance_and_magnetic_flux_density = calculate_inductance_and_magnetic_flux_density(magnetic, operatingPoint);

    return inductance_and_magnetic_flux_density.first;
}

int MagnetizingInductance::calculate_number_turns_from_gapping_and_inductance(Core core, Coil coil, Inputs* inputs, DimensionalValues preferredValue) {
    // Single source of truth for "turns for a target inductance": the EXACT
    // inverse of calculate_inductance_from_number_turns_and_gapping (the model
    // the inductance filter uses), so the two can never disagree. The previous
    // implementation ran its own DC-bias permeability-refinement loop that was
    // ill-conditioned for gapped cores under heavy bias — it drove the looked-up
    // permeability toward (and below) zero and over-counted N ~100x (a gapped
    // ferrite flyback E-core: N 13 -> 1631). Instead: seed from the unbiased
    // reluctance, then Newton-step on L proportional to N^2 against the canonical
    // operating-point inductance, then ensure it clears the (rounded) target.
    double desiredMagnetizingInductance = resolve_dimensional_values(
        inputs->get_design_requirements().get_magnetizing_inductance(), preferredValue);

    if (desiredMagnetizingInductance <= 0 || coil.get_functional_description().empty()) {
        // Lm = 0 means "not specified" (see pre_process_inputs); nothing to size.
        return std::max<int>(1, coil.get_functional_description().empty()
                                ? 1 : static_cast<int>(coil.get_functional_description()[0].get_number_turns()));
    }

    double frequency = Defaults().coreAdviserFrequencyReference;
    double temperature = Defaults().ambientTemperature;
    OperatingPoint operatingPoint;
    OperatingPoint* operatingPointPtr = nullptr;
    if (inputs->get_operating_points().size() > 0) {
        operatingPoint = inputs->get_operating_point(0);
        temperature = operatingPoint.get_conditions().get_ambient_temperature();
        frequency = operatingPoint.get_mutable_excitations_per_winding()[0].get_frequency();
        OperatingPointExcitation excitation = Inputs::get_primary_excitation(operatingPoint);
        if (!excitation.get_magnetizing_current() && !excitation.get_voltage()) {
            Inputs::set_current_as_magnetizing_current(&operatingPoint);
        }
        operatingPointPtr = &operatingPoint;
    }

    // Unbiased reluctance seed: N0 = round(sqrt(L * R_core(mu_initial))).
    OpenMagnetics::InitialPermeability initialPermeability;
    ReluctanceModels reluctanceModelEnum;
    from_json(_models["gapReluctance"], reluctanceModelEnum);
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(reluctanceModelEnum);
    double initialMu = initialPermeability.get_initial_permeability(core.resolve_material(), temperature, std::nullopt, frequency);
    double seedReluctance = reluctanceModel->get_core_reluctance(core, initialMu).get_core_reluctance();
    int numberTurnsPrimary = std::max(1, static_cast<int>(std::round(std::sqrt(desiredMagnetizingInductance * seedReluctance))));

    if (operatingPointPtr == nullptr) {
        // No operating point: the unbiased seed is the exact answer.
        return numberTurnsPrimary;
    }

    // Evaluate the canonical operating-point inductance at a trial turn count.
    auto inductanceAtTurns = [&](int n) -> double {
        coil.get_mutable_functional_description()[0].set_number_turns(static_cast<int64_t>(n));
        auto out = calculate_inductance_from_number_turns_and_gapping(core, coil, operatingPointPtr);
        return resolve_dimensional_values(out.get_magnetizing_inductance());
    };

    // Newton on L proportional to N^2: corrects the seed up or down.
    for (int iteration = 0; iteration < 6; ++iteration) {
        double inductance = inductanceAtTurns(numberTurnsPrimary);
        if (inductance <= 0) break;
        int next = std::max(1, static_cast<int>(std::round(numberTurnsPrimary * std::sqrt(desiredMagnetizingInductance / inductance))));
        if (next == numberTurnsPrimary) break;
        numberTurnsPrimary = next;
    }
    if (preferredValue == DimensionalValues::MINIMUM) {
        // The caller asked for a lower bound (e.g. the core adviser sizes
        // against the minimum requirement): the target is a hard floor, so
        // bump until the operating inductance actually clears it (integer
        // rounding and real permeability rolloff can leave it just under).
        for (int bump = 0; bump < 100; ++bump) {
            double inductance = inductanceAtTurns(numberTurnsPrimary);
            if (inductance <= 0 || inductance >= desiredMagnetizingInductance) break;
            numberTurnsPrimary += 1;
        }
    }
    else {
        // A nominal/typical target is not a hard floor: integer turns cannot
        // hit it exactly, so pick the neighbour with the smallest absolute
        // error. Ceil-style bumping here accepted a +26.6% overshoot to avoid
        // a -0.6% undershoot (ABT #600).
        double bestError = std::numeric_limits<double>::infinity();
        int bestTurns = numberTurnsPrimary;
        for (int candidate : {numberTurnsPrimary - 1, numberTurnsPrimary, numberTurnsPrimary + 1}) {
            if (candidate < 1) {
                continue;
            }
            double inductance = inductanceAtTurns(candidate);
            if (inductance <= 0) {
                continue;
            }
            double error = std::abs(inductance - desiredMagnetizingInductance);
            if (error < bestError) {
                bestError = error;
                bestTurns = candidate;
            }
        }
        numberTurnsPrimary = bestTurns;
    }

    return std::max(1, numberTurnsPrimary);
}

int MagnetizingInductance::calculate_number_turns_from_gapping_and_inductance(Core core, Inputs* inputs, DimensionalValues preferredValue) {
    // Legacy 3-argument shim (no coil) — see header. Builds the single-primary-
    // winding coil the old coil-less implementation implicitly assumed (one
    // winding on this core's bobbin) and forwards to the canonical coil-aware
    // overload. Magnetizing inductance referred to the primary (N^2 / R_core) is
    // independent of wire gauge and winding layout, so the placeholder winding
    // does not bias the result; the canonical version overwrites its turn count
    // while solving. Multi-winding callers must use the coil-taking overload.
    Coil coil;
    coil.set_bobbin(Bobbin::create_quick_bobbin(core));
    Winding winding;
    winding.set_name("winding 0");
    winding.set_number_turns(1);
    winding.set_number_parallels(1);
    winding.set_isolation_side_from_index(0);
    winding.set_wire("Round 0.475 - Grade 1");
    coil.get_mutable_functional_description().push_back(winding);
    // No wind() needed: the canonical path derives Lm from the core reluctance and
    // functional_description turn count, never from wound turn positions (the
    // coil-aware tests pass an unwound coil too).
    return calculate_number_turns_from_gapping_and_inductance(core, coil, inputs, preferredValue);
}

double MagnetizingInductance::calculate_gap_from_saturation_constraint(Core core,
                                                                       Inputs* inputs,
                                                                       double targetMagneticFluxDensity,
                                                                       double magnetizingCurrentPeak) {
    double desiredMagnetizingInductance = resolve_dimensional_values(inputs->get_design_requirements().get_magnetizing_inductance(), DimensionalValues::NOMINAL);
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    
    ReluctanceModels reluctanceModelEnum;
    from_json(_models["gapReluctance"], reluctanceModelEnum);
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(reluctanceModelEnum);
    
    // Start with energy-based gap as initial guess:
    // 0.5*L*I^2 = B^2/(2*mu0) * A * g  =>  g = mu0*L*I^2 / (A*B^2)
    // (a spurious factor 2 here used to double the initial guess, and since the
    // search below only ever GROWS the gap, it was returned as-is)
    auto constants = OpenMagnetics::Constants();
    double initialGap = (desiredMagnetizingInductance * pow(magnetizingCurrentPeak, 2) * constants.vacuumPermeability) /
                        (effectiveArea * pow(targetMagneticFluxDensity, 2));
    
    double gapLength = initialGap;
    const double gapStepFactor = 1.5;  // Increase gap by 50% per iteration
    const int maxIterations = 15;
    
    for (int iteration = 0; iteration < maxIterations; ++iteration) {
        // Create core with current gap
        Core testCore = core;
        auto columns = testCore.get_processed_description().value().get_columns();
        
        auto basicCentralGap = CoreGap();
        basicCentralGap.set_type(GapType::SUBTRACTIVE);
        basicCentralGap.set_length(gapLength);
        basicCentralGap.set_area(columns[0].get_area());
        
        auto basicLateralGap = CoreGap();
        basicLateralGap.set_type(GapType::RESIDUAL);
        basicLateralGap.set_length(constants.residualGap);
        
        std::vector<CoreGap> gapping;
        gapping.push_back(basicCentralGap);
        for (size_t i = 1; i < columns.size(); ++i) {
            auto lateralGap = basicLateralGap;
            lateralGap.set_area(columns[i].get_area());
            gapping.push_back(lateralGap);
        }
        testCore.get_mutable_functional_description().set_gapping(gapping);
        testCore.process_gap_or_throw();
        
        // Calculate total reluctance manually (classic formula)
        // Core reluctance: R_core = l_e / (μ₀ * μ_r * A_e), using the MATERIAL's
        // initial permeability (a hardcoded 2000 skewed the search for powder
        // cores with μr 14-125 and high-permeability MnZn ferrites alike)
        double effectiveLength = testCore.get_processed_description()->get_effective_parameters().get_effective_length();
        double initialPermeabilityValue = InitialPermeability::get_initial_permeability(testCore.resolve_material());
        double coreReluctance = effectiveLength / (constants.vacuumPermeability * initialPermeabilityValue * effectiveArea);
        
        // Gap reluctance: R_gap = l_gap / (μ₀ * A_gap)
        double gapArea = columns[0].get_area();
        double gapReluctance = gapLength / (constants.vacuumPermeability * gapArea);
        
        double totalReluctance = coreReluctance + gapReluctance;
        
        // Calculate turns for this gap to achieve target inductance
        int numberTurns = std::round(sqrt(desiredMagnetizingInductance * totalReluctance));
        numberTurns = std::max(1, numberTurns);
        
        // Calculate B-field: B = (N * I) / (R_total * A_e)
        double magneticFlux = numberTurns * magnetizingCurrentPeak / totalReluctance;
        double bPeak = magneticFlux / effectiveArea;
        
        // Check convergence
        if (bPeak <= targetMagneticFluxDensity) {
            return gapLength;
        }
        
        // Increase gap for next iteration
        gapLength *= gapStepFactor;
    }
    
    // Return last calculated gap even if not converged
    return gapLength;
}

Core get_core_with_ground_gapping(Core core, double gapLength) {
    auto constants = OpenMagnetics::Constants();
    auto basicCentralGap = CoreGap();
    basicCentralGap.set_type(GapType::SUBTRACTIVE);
    basicCentralGap.set_length(gapLength);
    auto basicLateralGap = CoreGap();
    basicLateralGap.set_type(GapType::RESIDUAL);
    basicLateralGap.set_length(constants.residualGap);
    std::vector<CoreGap> gapping;
    gapping.push_back(basicCentralGap);
    for (size_t i = 0; i < core.get_processed_description().value().get_columns().size() - 1; ++i) {
        gapping.push_back(basicLateralGap);
    }
    core.get_mutable_functional_description().set_gapping(gapping);
    core.process_gap_or_throw();
    return core;
}

Core get_core_with_distributed_gapping(Core core, double gapLength, size_t numberDistributedGaps) {
    auto constants = OpenMagnetics::Constants();
    auto basicCentralGap = CoreGap();
    basicCentralGap.set_type(GapType::SUBTRACTIVE);
    basicCentralGap.set_length(gapLength);
    auto basicLateralGap = CoreGap();
    basicLateralGap.set_type(GapType::RESIDUAL);
    basicLateralGap.set_length(constants.residualGap);
    std::vector<CoreGap> gapping;
    for (size_t i = 0; i < numberDistributedGaps; ++i) {
        gapping.push_back(basicCentralGap);
    }
    for (size_t i = 0; i < core.get_processed_description().value().get_columns().size() - 1; ++i) {
        gapping.push_back(basicLateralGap);
    }
    core.get_mutable_functional_description().set_gapping(gapping);
    core.process_gap_or_throw();
    return core;
}

Core get_core_with_spacer_gapping(Core core, double gapLength) {
    auto basicCentralGap = CoreGap();
    basicCentralGap.set_type(GapType::ADDITIVE);
    basicCentralGap.set_length(gapLength);
    auto basicLateralGap = CoreGap();
    basicLateralGap.set_type(GapType::ADDITIVE);
    basicLateralGap.set_length(gapLength);
    std::vector<CoreGap> gapping;
    gapping.push_back(basicCentralGap);
    for (size_t i = 0; i < core.get_processed_description().value().get_columns().size() - 1; ++i) {
        gapping.push_back(basicLateralGap);
    }
    core.get_mutable_functional_description().set_gapping(gapping);
    core.process_gap_or_throw();
    return core;
}

std::vector<CoreGap> MagnetizingInductance::calculate_gapping_from_number_turns_and_inductance(Core core,
                                                                                               Coil coil,
                                                                                               Inputs* inputs,
                                                                                               GappingType gappingType,
                                                                                               size_t decimals) {
    double frequency = Defaults().coreAdviserFrequencyReference;
    double temperature = Defaults().ambientTemperature;
    OperatingPointExcitation excitation;
    OperatingPoint operatingPoint;

    auto constants = OpenMagnetics::Constants();
    if (inputs->get_operating_points().size() > 0) {
        operatingPoint = inputs->get_operating_point(0);
        excitation = Inputs::get_primary_excitation(operatingPoint);
        temperature = operatingPoint.get_conditions().get_ambient_temperature();
        frequency = operatingPoint.get_mutable_excitations_per_winding()[0].get_frequency();
    }

    double numberTurnsPrimary = coil.get_functional_description()[0].get_number_turns();
    double desiredMagnetizingInductance = resolve_dimensional_values(inputs->get_design_requirements().get_magnetizing_inductance());

    // ABT #635: a caller may hand us a Core built from a bare functionalDescription (no
    // processedDescription yet) -- that is a legal MAS core and every neighbouring entry point
    // accepts it. Dereferencing the optional unconditionally turned that into an opaque
    // "bad optional access" for the whole API. Process it here instead, exactly as
    // calculate_core_maximum_magnetic_energy() and friends already do. `core` is taken BY VALUE,
    // so processing it is local to this call and cannot surprise the caller.
    if (!core.get_processed_description()) {
        core.process_data();
        core.process_gap_or_throw();
    }
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    OpenMagnetics::InitialPermeability initialPermeability;
    size_t timeout;
    double currentInitialPermeability;

    ReluctanceModels reluctanceModelEnum;
    from_json(_models["gapReluctance"], reluctanceModelEnum);
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(reluctanceModelEnum);
    double neededTotalReluctance = pow(numberTurnsPrimary, 2) / desiredMagnetizingInductance;

    currentInitialPermeability = initialPermeability.get_initial_permeability(core.resolve_material(), temperature, std::nullopt, frequency);

    if (!excitation.get_voltage() && excitation.get_current()) {
        Inputs::set_current_as_magnetizing_current(&operatingPoint);
        inputs->set_operating_point_by_index(operatingPoint, 0);
    }

    // DC bias (ABT #1093). The flux density the target inductance imposes is known
    // (B = N·I/(R·Ae) with the NEEDED reluctance: at the solution the gapped core has
    // exactly that reluctance), and the material permeability under that bias is its
    // reversible permeability at the field strength that carries B on the material's
    // own magnetisation curve. The former fixed-point iteration µ ← µ(B/(µ0·µ)) is a
    // secant through the incremental curve; it diverges wherever |dµ/dH|·B/(µ0·µ²) > 1,
    // the knee of any ferrite, ran away to a saturated permeability, and the search then
    // clamped at the residual gap for any target above a few tens of µH.
    if (excitation.get_magnetizing_current()) {
        auto magneticFlux = OpenMagnetics::MagneticField::calculate_magnetic_flux(operatingPoint.get_mutable_excitations_per_winding()[0].get_magnetizing_current().value(), neededTotalReluctance, numberTurnsPrimary);
        auto magneticFluxDensity = OpenMagnetics::MagneticField::calculate_magnetic_flux_density(magneticFlux, effectiveArea);
        double biasFluxDensity = fabs(magneticFluxDensity.get_processed().value().get_offset());
        if (biasFluxDensity > 0) {
            auto material = core.resolve_material();
            double biasFieldStrength = InitialPermeability::get_magnetic_field_dc_bias_for_flux_density(material, biasFluxDensity, temperature, frequency);
            currentInitialPermeability = initialPermeability.get_initial_permeability(material, temperature, biasFieldStrength, frequency);
        }
    }

    if (!excitation.get_voltage() && excitation.get_current()) {
        operatingPoint.get_mutable_excitations_per_winding()[0].set_voltage(Inputs::calculate_induced_voltage(excitation, desiredMagnetizingInductance));
        inputs->set_operating_point_by_index(operatingPoint, 0);
    }

    double gapLength = constants.residualGap;
    double gapLengthModification = constants.initialGapLengthForSearching;
    bool increasingGap = true;
    double fringingFactorOneGap = 0;

    double reluctance = 0;
    size_t numberDistributedGaps = 3;
    timeout = 100;
    Core gappedCore;

    while (true) {
        reluctance = 0;
        switch (gappingType) {
            case GappingType::GROUND:
                gappedCore = get_core_with_ground_gapping(core, gapLength);
                break;
            case GappingType::SPACER:
                gappedCore = get_core_with_spacer_gapping(core, gapLength);
                break;
            case GappingType::RESIDUAL:
                throw GapException("Residual type cannot be chosen to calculate the needed gapping");
                break;
            case GappingType::DISTRIBUTED:
                while (numberDistributedGaps > 3) {
                    gappedCore = get_core_with_distributed_gapping(core, gapLength, numberDistributedGaps);
                    fringingFactorOneGap = reluctanceModel->get_gap_reluctance(gappedCore.get_gapping()[0]).get_fringing_factor();
                    if (fringingFactorOneGap < constants.minimumDistributedFringingFactor &&
                        numberDistributedGaps > 1) {
                        gapLength *= numberDistributedGaps;
                        numberDistributedGaps -= 2;
                        gapLength /= numberDistributedGaps;
                    }
                    else {
                        break;
                    }
                }
                while (true) {
                    gappedCore = get_core_with_distributed_gapping(core, gapLength, numberDistributedGaps);
                    fringingFactorOneGap = reluctanceModel->get_gap_reluctance(gappedCore.get_gapping()[0]).get_fringing_factor();
                    if (fringingFactorOneGap > constants.maximumDistributedFringingFactor) {
                        gapLength *= numberDistributedGaps;
                        numberDistributedGaps += 2;
                        gapLength /= numberDistributedGaps;
                    }
                    else {
                        break;
                    }
                }
                break;
            default:
                throw GapException("Unknown type of gap, options are {GROUND, SPACER, RESIDUAL, DISTRIBUTED}");
        }

        auto magnetizingInductanceOutput = reluctanceModel->get_core_reluctance(gappedCore, currentInitialPermeability);
        reluctance = magnetizingInductanceOutput.get_core_reluctance();

        if (fabs(neededTotalReluctance - reluctance) / neededTotalReluctance < 0.001 || timeout == 0) {
            break;
        }
        else {
            if (neededTotalReluctance < reluctance && increasingGap) {
                increasingGap = false;
                gapLengthModification = std::max(gapLengthModification / 2., constants.residualGap);
            }
            if (neededTotalReluctance > reluctance && !increasingGap) {
                increasingGap = true;
                gapLengthModification = std::max(gapLengthModification / 2., constants.residualGap);
            }
            if (increasingGap) {
                gapLength += gapLengthModification;
            }
            else {
                gapLength = std::max(constants.residualGap, gapLength - gapLengthModification);
            }

            timeout--;
        }
    }

    gapLength = std::max(constants.residualGap, roundFloat(gapLength, decimals));

    switch (gappingType) {
        case GappingType::GROUND:
            return get_core_with_ground_gapping(core, gapLength).get_gapping();
        case GappingType::SPACER:
            return get_core_with_spacer_gapping(core, gapLength).get_gapping();
        case GappingType::RESIDUAL:
            throw GapException("Residual type cannot be chosen to calculate the needed gapping");
            break;
        case GappingType::DISTRIBUTED:
            return get_core_with_distributed_gapping(core, gapLength, numberDistributedGaps).get_gapping();
            break;
        default:
            throw GapException("Unknown type of gap, options are {GROUND, SPACER, RESIDUAL, DISTRIBUTED}");
    }
}

double MagnetizingInductance::calculate_turns_for_gap(
    Core& core,
    double targetInductance,
    double temperature,
    double frequency) 
{
    OpenMagnetics::InitialPermeability initialPermeability;
    ReluctanceModels reluctanceModelEnum;
    from_json(_models["gapReluctance"], reluctanceModelEnum);
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(reluctanceModelEnum);
    
    double currentInitialPermeability = initialPermeability.get_initial_permeability(
        core.resolve_material(), temperature, std::nullopt, frequency);
    
    auto magnetizingInductanceOutput = reluctanceModel->get_core_reluctance(core, currentInitialPermeability);
    double totalReluctance = magnetizingInductanceOutput.get_core_reluctance();
    
    double turns = sqrt(targetInductance * totalReluctance);
    return std::max(1.0, std::round(turns));
}

double MagnetizingInductance::calculate_flux_density_peak(
    Core& core,
    double numberTurns,
    double peakCurrent,
    double temperature,
    double frequency) 
{
    OpenMagnetics::InitialPermeability initialPermeability;
    ReluctanceModels reluctanceModelEnum;
    from_json(_models["gapReluctance"], reluctanceModelEnum);
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(reluctanceModelEnum);
    
    double currentInitialPermeability = initialPermeability.get_initial_permeability(
        core.resolve_material(), temperature, std::nullopt, frequency);
    
    auto magnetizingInductanceOutput = reluctanceModel->get_core_reluctance(core, currentInitialPermeability);
    double totalReluctance = magnetizingInductanceOutput.get_core_reluctance();
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    
    // B = Φ / A = (N * I / R) / A
    double magneticFlux = numberTurns * peakCurrent / totalReluctance;
    double bPeak = magneticFlux / effectiveArea;
    
    return bPeak;
}

double MagnetizingInductance::calculate_flux_density_peak_from_voltage(
    Core& core,
    double numberTurns,
    double voltagePeak,
    double frequency) 
{
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    
    // Faraday's Law: V = N * Ae * dB/dt
    // For sinusoidal: V_peak = N * Ae * ω * B_peak
    // Therefore: B_peak = V_peak / (N * Ae * ω)
    double omega = 2 * std::numbers::pi * frequency;
    
    if (numberTurns <= 0 || effectiveArea <= 0 || omega <= 0) {
        return 0.0;
    }
    
    double bPeak = voltagePeak / (numberTurns * effectiveArea * omega);
    return bPeak;
}

double MagnetizingInductance::calculate_flux_density_peak_from_volt_seconds(
    Core& core,
    double numberTurns,
    double maxVoltSeconds)
{
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    if (numberTurns <= 0 || effectiveArea <= 0 || maxVoltSeconds <= 0) {
        return 0.0;
    }
    // Faraday's Law: V = N · A_e · dB/dt  →  B_peak = max|∫V dt| / (N · A_e).
    // Works for arbitrary waveforms (DAB square, CLLLC trapezoid, ...);
    // does NOT assume sinusoidal.
    return maxVoltSeconds / (numberTurns * effectiveArea);
}

double MagnetizingInductance::calculate_turns_from_volt_seconds_and_max_flux_density(
    Core& core,
    double maxVoltSeconds,
    double maxFluxDensity)
{
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    if (effectiveArea <= 0 || maxFluxDensity <= 0 || maxVoltSeconds <= 0) {
        return 0.0;
    }
    // Inverse of B_peak = V·s / (N · A_e): N_min = ceil(V·s / (B · A_e)).
    return std::max(1.0, std::ceil(maxVoltSeconds / (maxFluxDensity * effectiveArea)));
}

double MagnetizingInductance::calculate_turns_from_voltage_and_max_flux_density(
    Core& core,
    double voltagePeak,
    double frequency,
    double maxFluxDensity)
{
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    
    // From Faraday's Law: B_peak = V_peak / (N * Ae * ω)
    // Therefore: N = V_peak / (Ae * ω * B_max)
    double omega = 2 * std::numbers::pi * frequency;
    
    if (effectiveArea <= 0 || omega <= 0 || maxFluxDensity <= 0) {
        return 1.0;
    }
    
    double turns = voltagePeak / (effectiveArea * omega * maxFluxDensity);
    return std::max(1.0, std::ceil(turns));
}

std::pair<double, double> MagnetizingInductance::calculate_optimal_gap_and_turns(
    const Core& core,
    Inputs* inputs,
    double maxFluxDensity,
    double peakCurrent) 
{
    auto constants = OpenMagnetics::Constants();
    double targetInductance = resolve_dimensional_values(
        inputs->get_design_requirements().get_magnetizing_inductance(), 
        DimensionalValues::NOMINAL);
    
    double temperature = Defaults().ambientTemperature;
    double frequency = Defaults().coreAdviserFrequencyReference;
    if (inputs->get_operating_points().size() > 0) {
        temperature = inputs->get_operating_point(0).get_conditions().get_ambient_temperature();
        frequency = inputs->get_operating_point(0).get_excitations_per_winding()[0].get_frequency();
    }
    
    // Create a mutable copy for processing
    Core tempCore = core;
    if (!tempCore.get_processed_description()) {
        tempCore.process_data();
    }
    
    double effectiveArea = tempCore.get_processed_description()->get_effective_parameters().get_effective_area();
    
    // Step 1: Calculate minimum turns from saturation constraint
    // B_max = L * I_peak / (N * A)  =>  N_min = L * I_peak / (B_max * A)
    double turnsFromSaturation = (targetInductance * peakCurrent) / (maxFluxDensity * effectiveArea);
    double minTurns = std::max(1.0, std::ceil(turnsFromSaturation));
    
    // Step 2: Calculate required reluctance for these turns
    // L = N² / R  =>  R = N² / L
    double requiredReluctance = (minTurns * minTurns) / targetInductance;
    
    // Step 3: Get core reluctance (ungapped) to determine needed gap reluctance
    OpenMagnetics::InitialPermeability initialPermeability;
    ReluctanceModels reluctanceModelEnum;
    from_json(_models["gapReluctance"], reluctanceModelEnum);
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(reluctanceModelEnum);
    
    double currentInitialPermeability = initialPermeability.get_initial_permeability(
        tempCore.resolve_material(), temperature, std::nullopt, frequency);
    
    // Estimate core reluctance from effective parameters
    double effectiveLength = tempCore.get_processed_description()->get_effective_parameters().get_effective_length();
    double coreReluctance = effectiveLength / (constants.vacuumPermeability * currentInitialPermeability * effectiveArea);
    
    // Step 4: Calculate needed gap reluctance
    double neededGapReluctance = requiredReluctance - coreReluctance;
    
    if (neededGapReluctance <= 0) {
        // No gap needed - core reluctance alone is enough
        // But verify saturation constraint is met
        double actualTurns = sqrt(targetInductance * coreReluctance);
        return {0.0, std::max(1.0, std::round(actualTurns))};
    }
    
    // Step 5: Calculate gap length from gap reluctance
    // R_gap = g / (μ₀ * A_gap)  =>  g = R_gap * μ₀ * A_gap
    // Note: This is simplified - doesn't account for fringing
    double gapArea = tempCore.get_processed_description()->get_columns()[0].get_area();
    double gapLength = neededGapReluctance * constants.vacuumPermeability * gapArea;
    
    // Step 6: Validate gap is practical
    double columnWidth = tempCore.get_columns()[0].get_width();
    double maxGap = 0.4 * columnWidth;  // 40% of column width max
    double minGap = constants.residualGap;
    
    if (gapLength > maxGap) {
        // Gap too large - core is too small for this application
        return {-1.0, -1.0};
    }
    
    gapLength = std::max(minGap, gapLength);
    
    // Step 7: Refine turns with actual gap (accounts for fringing via reluctance model)
    Core gappedCore = tempCore;
    gappedCore.set_ground_gapping(gapLength);
    gappedCore.process_gap_or_throw();
    
    double finalTurns = calculate_turns_for_gap(gappedCore, targetInductance, temperature, frequency);
    
    // Step 8: Verify saturation constraint with final values
    double finalBpeak = calculate_flux_density_peak(gappedCore, finalTurns, peakCurrent, temperature, frequency);
    
    if (finalBpeak > maxFluxDensity * 1.05) {  // Allow 5% tolerance
        // Need to increase turns to reduce B
        finalTurns = std::ceil((targetInductance * peakCurrent) / (maxFluxDensity * effectiveArea));
        // Recalculate gap for new turns
        double newRequiredReluctance = (finalTurns * finalTurns) / targetInductance;
        double newNeededGapReluctance = newRequiredReluctance - coreReluctance;
        if (newNeededGapReluctance > 0) {
            gapLength = newNeededGapReluctance * constants.vacuumPermeability * gapArea;
            gapLength = std::min(maxGap, std::max(minGap, gapLength));
        }
    }
    
    return {gapLength, finalTurns};
}

} // namespace OpenMagnetics
