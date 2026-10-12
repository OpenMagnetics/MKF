#pragma once
#include "Models.h"
#include "constructive_models/Magnetic.h"
#include <map>
#include <optional>
#include <string>
#include <vector>

using namespace MAS;

namespace OpenMagnetics {

// ABT #1304: the cores a winding-loss model set can model. The field models are not interchangeable across these:
// the image lattice needs one rectangular window, the Kelvin images a toroid, ALBACH a round central post.
enum class WindingLossesModelSetScope : int {
    ANY,                        // every core
    SINGLE_RECTANGULAR_WINDOW,  // one winding window, not a toroid
    ROUND_POST_SINGLE_WINDOW,   // one winding window around a round central column (PQ, RM, P, EP, ETD, ...)
    TOROID,
    MULTI_WINDOW,               // several winding windows (multi-column cores)
    NO_WINDOW_FRAME             // drum / rod / molded families (DRUM, ROD, DRUM_RING, DRUM_SEMISHIELDED, DRUM_PLATE, MOLDED):
                                // no ferrite frame closes an air window, so the neutral-window image lattice has nothing to stand on
};

// A named, self-consistent combination of the winding-loss models: ONE field solution for the whole magnetic (field
// model + fringing/core model) and, per wire type, the skin and proximity models that consume it. Models that only
// make sense together are chosen together, and a set refuses -- never silently adapts to -- a magnetic it was not
// built for.
struct WindingLossesModelSet {
    WindingLossesModelSets id;
    std::string name;         // human-readable, as UIs show it
    std::string description;  // what it is and when to use it
    std::string validation;   // measured evidence and its reference, or "not validated"
    // nullopt: MKF's own per-wire-type choice of field model (the historical automatic rule, LEGACY_DEFAULT only).
    std::optional<MagneticFieldStrengthModels> fieldModel;
    MagneticFieldStrengthFringingEffectModels fringingModel;
    std::map<WireType, WindingSkinEffectLossesModels> skinModelPerWireType;
    std::map<WireType, WindingProximityEffectLossesModels> proximityModelPerWireType;
    std::vector<WindingLossesModelSetScope> scopes;
    bool forDesign;           // false: benchmark / legacy baseline, never a design reference
};

// The registry. CUSTOM is not in it: CUSTOM means "the individually selected models" (Settings / WindingLossesModels).
const WindingLossesModelSet& get_winding_losses_model_set(WindingLossesModelSets id);
std::vector<WindingLossesModelSets> get_all_winding_losses_model_sets();

// Throws ModelNotAvailableException naming what does not fit: the core (scope) or a wire type without models in the set.
void check_winding_losses_model_set_applies(const WindingLossesModelSet& set, Magnetic magnetic);
// The sets that apply to this magnetic (CUSTOM always does).
std::vector<WindingLossesModelSets> get_applicable_winding_losses_model_sets(Magnetic magnetic);

// MKF's historical automatic choice of field model by wire type: the CUSTOM path and LEGACY_DEFAULT (ABT #1304).
MagneticFieldStrengthModels automatic_winding_losses_field_model(Magnetic& magnetic);
// The field model a set solves with: its own, or the historical automatic choice when it names none (LEGACY_DEFAULT).
MagneticFieldStrengthModels resolve_winding_losses_field_model(const WindingLossesModelSet& set, Magnetic& magnetic);

// ABT #1759: the set used when no winding-loss model is chosen. Per core type, the set most accurate against 2D OMFEM
// (Defaults.h): the neutral window where its frame exists (one winding window, not a toroid, every turn crossing the
// section once), the historical automatic choice elsewhere.
WindingLossesModelSets get_default_winding_losses_model_set(Magnetic magnetic);
// The set Settings select for this magnetic, for every consumer of the winding-loss field (winding losses, the Painter's
// field plot). Settings' set when one is named. Otherwise (CUSTOM) the default set for the core while the field
// configuration is the default one -- field and fringing models at their Defaults and fringing included -- and CUSTOM
// (the individually configured models) once any of them is configured.
WindingLossesModelSets resolve_winding_losses_model_set_from_settings(Magnetic magnetic, bool includeFringing);

} // namespace OpenMagnetics
