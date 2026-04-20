#pragma once
#include <cstddef>
#include <utility>
#include <vector>

namespace OpenMagnetics {

// Returns the precomputed CCI strand coordinates for N strands, normalized to [-1, 1].
// Returns nullptr if N > 1000 (caller should fall back to spiral algorithm).
const std::vector<std::pair<float,float>>* get_cci_coordinates(size_t n);

} // namespace OpenMagnetics
