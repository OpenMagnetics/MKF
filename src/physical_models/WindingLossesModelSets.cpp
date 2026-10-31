#include "physical_models/WindingLossesModelSets.h"
#include "support/Exceptions.h"
#include "support/Settings.h"
#include "Defaults.h"
#include <magic_enum.hpp>

namespace OpenMagnetics {

namespace {

using Field = MagneticFieldStrengthModels;
using Fringing = MagneticFieldStrengthFringingEffectModels;
using Skin = WindingSkinEffectLossesModels;
using Proximity = WindingProximityEffectLossesModels;
using Scope = WindingLossesModelSetScope;

// MKF's historical per-wire skin choice (WindingSkinEffectLosses::get_model), kept identical by every set that does
// not deliberately change it.
const std::map<WireType, Skin> historicalSkin = {
    {WireType::ROUND, Skin::ALBACH}, {WireType::LITZ, Skin::ALBACH}, {WireType::PLANAR, Skin::ALBACH},
    {WireType::RECTANGULAR, Skin::KUTKUT}, {WireType::FOIL, Skin::KUTKUT}};

std::map<WindingLossesModelSets, WindingLossesModelSet> build_registry() {
    std::map<WindingLossesModelSets, WindingLossesModelSet> registry;
    auto add = [&](WindingLossesModelSet set) { registry.emplace(set.id, std::move(set)); };

    add({WindingLossesModelSets::LEGACY_DEFAULT, "Legacy default",
         "MKF's historical automatic choice: field model by wire type (LAMMERANER round/litz, BINNS_LAWRENSON "
         "rectangular, WANG planar and foil), ROSHEN gap fringing, FERREIRA proximity for round/litz, WANG for flat wire.",
         "Reference behaviour. Known limits: a flat turn's current is counted twice under BINNS_LAWRENSON "
         "with the Wang mesher, residual gaps carry no field, the image lattice holds net current.",
         std::nullopt, Fringing::ROSHEN, historicalSkin,
         {{WireType::ROUND, Proximity::FERREIRA}, {WireType::LITZ, Proximity::FERREIRA}, {WireType::PLANAR, Proximity::WANG},
          {WireType::RECTANGULAR, Proximity::WANG}, {WireType::FOIL, Proximity::WANG}},
         {Scope::ANY}, true});

    add({WindingLossesModelSets::NEUTRAL_WINDOW, "Neutral window (MMF boundary + Martinez)",
         "Charge-neutral window field (MMF_BOUNDARY on BINNS_LAWRENSON): the core's gap and ferrite MMF drops on the "
         "window boundary, every conductor imaged with its own images. MARTINEZ field-exclusion proximity for flat wire, "
         "FERREIRA (exact Bessel) for round and litz.",
         "Flat wire, per-turn R_ac/R_dc vs 2D OMFEM: E 42/21/20 0.75-1.73 per turn, 0.92-0.98 in total; PQ 27/17 "
         "0.67-1.29 / 0.79-0.82. Round and litz: not yet validated per turn.",
         Field::BINNS_LAWRENSON, Fringing::MMF_BOUNDARY, historicalSkin,
         {{WireType::ROUND, Proximity::FERREIRA}, {WireType::LITZ, Proximity::FERREIRA}, {WireType::PLANAR, Proximity::MARTINEZ},
          {WireType::RECTANGULAR, Proximity::MARTINEZ}, {WireType::FOIL, Proximity::MARTINEZ}},
         {Scope::SINGLE_RECTANGULAR_WINDOW}, true});

    add({WindingLossesModelSets::NEUTRAL_WINDOW_WANG, "Neutral window with Wang proximity",
         "The neutral-window field with WANG slab proximity for flat wire, for comparison with NEUTRAL_WINDOW.",
         "Flat wire vs 2D OMFEM: E 42/21/20 totals 0.78-0.89. Round and litz: not validated per turn.",
         Field::BINNS_LAWRENSON, Fringing::MMF_BOUNDARY, historicalSkin,
         {{WireType::ROUND, Proximity::FERREIRA}, {WireType::LITZ, Proximity::FERREIRA}, {WireType::PLANAR, Proximity::WANG},
          {WireType::RECTANGULAR, Proximity::WANG}, {WireType::FOIL, Proximity::WANG}},
         {Scope::SINGLE_RECTANGULAR_WINDOW}, true});

    add({WindingLossesModelSets::WANG_PLANAR, "Wang planar",
         "Wang's field model with its own two-point mesh and Wang proximity: the original calibrated pair for flat "
         "conductors.",
         "Field and loss are a calibrated pair: their product scales as width^2/thickness, right near the geometry it "
         "was tuned on.",
         Field::WANG, Fringing::ROSHEN,
         {{WireType::PLANAR, Skin::ALBACH}, {WireType::RECTANGULAR, Skin::KUTKUT}, {WireType::FOIL, Skin::KUTKUT}},
         {{WireType::PLANAR, Proximity::WANG}, {WireType::RECTANGULAR, Proximity::WANG}, {WireType::FOIL, Proximity::WANG}},
         {Scope::SINGLE_RECTANGULAR_WINDOW}, true});

    add({WindingLossesModelSets::SULLIVAN_GAP, "Sullivan gap",
         "LAMMERANER conductor field with SULLIVAN gap fringing (gap MMF as imaged filaments), FERREIRA proximity.",
         "Not validated per turn. Its gap sources are not neutral against the conductor lattice.",
         Field::LAMMERANER, Fringing::SULLIVAN, {{WireType::ROUND, Skin::ALBACH}, {WireType::LITZ, Skin::ALBACH}},
         {{WireType::ROUND, Proximity::FERREIRA}, {WireType::LITZ, Proximity::FERREIRA}},
         {Scope::SINGLE_RECTANGULAR_WINDOW}, true});

    add({WindingLossesModelSets::AXISYMMETRIC, "Axisymmetric (Albach)",
         "ALBACH axisymmetric conductor field around a round central post with ALBACH gap fringing (gaps outside its "
         "fitted range go through Roshen), FERREIRA proximity.",
         "Not validated per turn. Albach's fitted equivalent current loop diverged on real gaps, so every "
         "gap in its validity range is currently routed to Roshen.",
         Field::ALBACH, Fringing::ALBACH, {{WireType::ROUND, Skin::ALBACH}, {WireType::LITZ, Skin::ALBACH}},
         {{WireType::ROUND, Proximity::FERREIRA}, {WireType::LITZ, Proximity::FERREIRA}},
         {Scope::ROUND_POST_SINGLE_WINDOW}, true});

    add({WindingLossesModelSets::TOROID, "Toroid",
         "BINNS_LAWRENSON conductor field with Kelvin images for the ring, ROSHEN fringing; FERREIRA proximity for "
         "round and litz, MARTINEZ for flat wire.",
         "Not validated per turn; flat wire on toroids never compared with FEM.",
         Field::BINNS_LAWRENSON, Fringing::ROSHEN, historicalSkin,
         {{WireType::ROUND, Proximity::FERREIRA}, {WireType::LITZ, Proximity::FERREIRA}, {WireType::PLANAR, Proximity::MARTINEZ},
          {WireType::RECTANGULAR, Proximity::MARTINEZ}, {WireType::FOIL, Proximity::MARTINEZ}},
         {Scope::TOROID}, true});

    add({WindingLossesModelSets::MULTI_WINDOW, "Multi-window",
         "BINNS_LAWRENSON conductor field with a mirror lattice per winding window, ROSHEN fringing; FERREIRA proximity "
         "for round and litz, MARTINEZ for flat wire.",
         "Not validated per turn. MMF_BOUNDARY is not extended to multi-window cores yet.",
         Field::BINNS_LAWRENSON, Fringing::ROSHEN, historicalSkin,
         {{WireType::ROUND, Proximity::FERREIRA}, {WireType::LITZ, Proximity::FERREIRA}, {WireType::PLANAR, Proximity::MARTINEZ},
          {WireType::RECTANGULAR, Proximity::MARTINEZ}, {WireType::FOIL, Proximity::MARTINEZ}},
         {Scope::MULTI_WINDOW}, true});

    add({WindingLossesModelSets::DOWELL_1D, "Dowell 1-D (baseline)",
         "Dowell's one-dimensional layer field with Dowell skin and proximity for every wire type.",
         "Legacy baseline only: a 1-D infinite-layer approximation. Never use it as a reference.",
         Field::DOWELL, Fringing::ROSHEN,
         {{WireType::ROUND, Skin::DOWELL}, {WireType::LITZ, Skin::DOWELL}, {WireType::PLANAR, Skin::DOWELL},
          {WireType::RECTANGULAR, Skin::DOWELL}, {WireType::FOIL, Skin::DOWELL}},
         {{WireType::ROUND, Proximity::DOWELL}, {WireType::LITZ, Proximity::DOWELL}, {WireType::PLANAR, Proximity::DOWELL},
          {WireType::RECTANGULAR, Proximity::DOWELL}, {WireType::FOIL, Proximity::DOWELL}},
         {Scope::SINGLE_RECTANGULAR_WINDOW}, false});

    add({WindingLossesModelSets::EWALD_BENCHMARK, "Ewald & Biela benchmark",
         "The neutral-window field with Ewald & Biela's homogeneous-surface-field proximity model for flat wire "
         "(EPE'23), FERREIRA for round and litz. For comparing against the published state of the art.",
         "Benchmark only. Against 2D FEM on single stadium conductors at 2:1 it spans 0.17-1.15 (median 0.44).",
         Field::BINNS_LAWRENSON, Fringing::MMF_BOUNDARY, historicalSkin,
         {{WireType::ROUND, Proximity::FERREIRA}, {WireType::LITZ, Proximity::FERREIRA}, {WireType::PLANAR, Proximity::EWALD},
          {WireType::RECTANGULAR, Proximity::EWALD}, {WireType::FOIL, Proximity::EWALD}},
         {Scope::SINGLE_RECTANGULAR_WINDOW}, false});

    return registry;
}

const std::map<WindingLossesModelSets, WindingLossesModelSet>& registry() {
    static const auto instance = build_registry();
    return instance;
}

Scope scope_of(Magnetic& magnetic) {
    auto core = magnetic.get_core();
    const auto family = core.get_shape_family();
    if (family == CoreShapeFamily::T) {
        return Scope::TOROID;
    }
    // ABT #1783 (Alf 2026-10-12): a drum, a rod, and a drum closed by a ring, a shell or a plate wind in a window that
    // no ferrite frame closes on all four sides (the drum-ring's annular clearance gap lies outside the window's height,
    // so it opens nowhere onto it). The neutral window places the core's MMF on that frame, so it has no frame to place it on.
    // A MOLDED core (Alf 2026-10-12) has its winding embedded in the magnetic material itself: there is no air window
    // bounded by a high-permeability frame either.
    if (family == CoreShapeFamily::DRUM || family == CoreShapeFamily::ROD || family == CoreShapeFamily::DRUM_RING ||
        family == CoreShapeFamily::DRUM_SEMISHIELDED || family == CoreShapeFamily::DRUM_PLATE ||
        family == CoreShapeFamily::MOLDED) {
        return Scope::NO_WINDOW_FRAME;
    }
    if (!core.get_processed_description()) {
        core.process_data();
    }
    if (core.get_winding_windows().size() > 1) {
        return Scope::MULTI_WINDOW;
    }
    auto columns = core.get_columns();
    if (!columns.empty() && columns[core.get_main_column_index()].get_shape() == ColumnShape::ROUND) {
        return Scope::ROUND_POST_SINGLE_WINDOW;
    }
    return Scope::SINGLE_RECTANGULAR_WINDOW;
}

bool scope_covers(Scope allowed, Scope actual) {
    if (allowed == Scope::ANY || allowed == actual) {
        return true;
    }
    // A round central post is still one rectangular window in cross-section.
    return allowed == Scope::SINGLE_RECTANGULAR_WINDOW && actual == Scope::ROUND_POST_SINGLE_WINDOW;
}

} // namespace

const WindingLossesModelSet& get_winding_losses_model_set(WindingLossesModelSets id) {
    if (id == WindingLossesModelSets::CUSTOM) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "CUSTOM is not a registered model set: it means the individually selected winding-loss models");
    }
    auto it = registry().find(id);
    if (it == registry().end()) {
        throw ModelNotAvailableException("Unknown winding-loss model set " + std::string(magic_enum::enum_name(id)));
    }
    return it->second;
}

std::vector<WindingLossesModelSets> get_all_winding_losses_model_sets() {
    std::vector<WindingLossesModelSets> all;
    for (auto id : magic_enum::enum_values<WindingLossesModelSets>()) {
        all.push_back(id);
    }
    return all;
}

void check_winding_losses_model_set_applies(const WindingLossesModelSet& set, Magnetic magnetic) {
    auto actual = scope_of(magnetic);
    bool covered = false;
    for (auto allowed : set.scopes) {
        covered = covered || scope_covers(allowed, actual);
    }
    if (!covered) {
        throw ModelNotAvailableException("Winding-loss model set " + std::string(magic_enum::enum_name(set.id)) +
                                         " is not available for this core family (" +
                                         std::string(magic_enum::enum_name(magnetic.get_core().get_shape_family())) +
                                         ", scope " + std::string(magic_enum::enum_name(actual)) + ")");
    }
    auto coil = magnetic.get_coil();
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        auto wireType = coil.get_wire_type(windingIndex);
        if (!set.skinModelPerWireType.contains(wireType) || !set.proximityModelPerWireType.contains(wireType)) {
            throw ModelNotAvailableException("Winding-loss model set " + std::string(magic_enum::enum_name(set.id)) +
                                             " has no models for " + std::string(magic_enum::enum_name(wireType)) +
                                             " wire (winding " + coil.get_functional_description()[windingIndex].get_name() + ")");
        }
    }
}

MagneticFieldStrengthModels automatic_winding_losses_field_model(Magnetic& magnetic) {
    bool isPlanar = true;
    bool isRectangular = true;
    bool isFoil = true;
    for (auto wire : magnetic.get_mutable_coil().get_wires()) {
        isPlanar = isPlanar && wire.get_type() == WireType::PLANAR;
        isRectangular = isRectangular && wire.get_type() == WireType::RECTANGULAR;
        isFoil = isFoil && wire.get_type() == WireType::FOIL;
    }
    if (isPlanar) {
        return Field::WANG;
    }
    if (isRectangular) {
        return Field::BINNS_LAWRENSON;
    }
    if (isFoil) {
        return Field::WANG;
    }
    for (auto wire : magnetic.get_mutable_coil().get_wires()) {
        if (wire.get_type() == WireType::FOIL) {
            // Mixed foil + round/litz (e.g. foil primary, litz secondary): the Wang mesher represents a foil inducing
            // turn as 2 end filaments, which only the Wang label-pair formulas can reconstruct into a sheet field.
            // Feeding them to LAMMERANER places the whole foil current at the foil ENDS and inverts the MMF profile
            // across the stack (ABT #182). ALBACH subdivides flat conductors into distributed filaments and handles
            // mixed conductor types.
            return Field::ALBACH;
        }
    }
    return Field::LAMMERANER;
}

MagneticFieldStrengthModels resolve_winding_losses_field_model(const WindingLossesModelSet& set, Magnetic& magnetic) {
    if (set.fieldModel) {
        return set.fieldModel.value();
    }
    return automatic_winding_losses_field_model(magnetic);
}

WindingLossesModelSets get_default_winding_losses_model_set(Magnetic magnetic) {
    const Defaults defaults;
    switch (scope_of(magnetic)) {
        case Scope::TOROID:
            return defaults.windingLossesModelSetDefaultToroid;
        case Scope::MULTI_WINDOW:
            return defaults.windingLossesModelSetDefaultMultiWindow;
        case Scope::NO_WINDOW_FRAME:
            return defaults.windingLossesModelSetDefaultNoWindowFrame;
        case Scope::SINGLE_RECTANGULAR_WINDOW:
        case Scope::ROUND_POST_SINGLE_WINDOW:
            break;
        case Scope::ANY:
            throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "A core's scope is never ANY");
    }
    // A turn with a second crossing (a winding on a lateral column) is not one conductor of the single window frame.
    const auto& turnsDescription = magnetic.get_coil().get_turns_description();
    if (turnsDescription) {
        const auto turns = turnsDescription.value();
        for (const auto& turn : turns) {
            if (turn.get_additional_coordinates()) {
                return defaults.windingLossesModelSetDefaultSecondCrossing;
            }
        }
    }
    return defaults.windingLossesModelSetDefaultSingleWindow;
}

WindingLossesModelSets resolve_winding_losses_model_set_from_settings(Magnetic magnetic, bool includeFringing) {
    auto& settings = Settings::GetInstance();
    if (settings.get_winding_losses_model_set() != WindingLossesModelSets::CUSTOM) {
        return settings.get_winding_losses_model_set();
    }
    // ABT #2169: skin and proximity models chosen in Settings are a choice like a field model is; a default set would
    // fix every model and ignore them.
    if (settings.get_coil_enable_user_winding_losses_models()) {
        return WindingLossesModelSets::CUSTOM;
    }
    const Defaults defaults;
    const bool defaultFieldConfiguration = includeFringing &&
        settings.get_magnetic_field_strength_model() == defaults.magneticFieldStrengthModelDefault &&
        settings.get_magnetic_field_strength_fringing_effect_model() == defaults.magneticFieldStrengthFringingEffectModelDefault;
    if (!defaultFieldConfiguration) {
        return WindingLossesModelSets::CUSTOM;
    }
    return get_default_winding_losses_model_set(magnetic);
}

std::vector<WindingLossesModelSets> get_applicable_winding_losses_model_sets(Magnetic magnetic) {
    std::vector<WindingLossesModelSets> applicable = {WindingLossesModelSets::CUSTOM};
    for (const auto& [id, set] : registry()) {
        try {
            check_winding_losses_model_set_applies(set, magnetic);
            applicable.push_back(id);
        }
        catch (const ModelNotAvailableException&) {
            // Not applicable is the answer being computed here, not an error to hide.
        }
    }
    return applicable;
}

} // namespace OpenMagnetics
