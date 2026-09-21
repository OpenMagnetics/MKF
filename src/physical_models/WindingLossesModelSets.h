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
    MULTI_WINDOW                // several winding windows (multi-column cores)
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

} // namespace OpenMagnetics
