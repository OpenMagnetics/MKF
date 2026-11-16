#pragma once

#include "constructive_models/Coil.h"

#include <optional>
#include <string>

namespace OpenMagnetics {

// ABT #2274 (Alf, 2026-11-05: "MKF counts the bump"): why a real-wound coil does not fit its core once
// its turns ride over the Z returns (Coil::judge_ride_over_core_contacts), or nothing when they fit. A
// core whose legs MKF's core room does not know, or a toroid with returns to ride over, cannot be shown
// to fit: its COIL_WINDING_ERROR comes back as the refusal. Shared by real_winding_refusal (the stored
// form) and the coil adviser's wind with planned leads, so both refuse with the same words.
std::optional<std::string> ride_over_core_refusal(Coil& coil);

} // namespace OpenMagnetics
