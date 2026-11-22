#pragma once

#include <MAS.hpp>
#include <vector>

namespace OpenMagnetics {

// ABT #2470 (Alf, 2026-11-22): the advanced materials file's B-H data enters MKF only when a caller loads that file
// (load_advanced_core_materials, as for its loss points), and in the layout the magnetisation curve reads: per
// temperature the rising branch up to the tip then the falling one back, from B = 0 and thinned (the file holds up to
// 8666 raw samples per temperature, which every material lookup would copy), or a single curve with H rising.
// A temperature no curve can be built from is left out, with a WARNING that says why.
std::vector<MAS::BhCycleElement> normalise_bh_cycle(const MAS::CoreMaterial& coreMaterial, const std::vector<MAS::BhCycleElement>& cycle);

} // namespace OpenMagnetics
