#pragma once
#include "MAS.hpp"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Wire.h"
#include "support/Utils.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <streambuf>
#include <vector>

using namespace MAS;

namespace OpenMagnetics {


enum class CoilMesherModels : int {
    WANG,
    CENTER
};


class CoilMesher {
  private:
  protected:
    double _quickModeForManyHarmonicsThreshold = 1;
  public:
    std::vector<Field> generate_mesh_inducing_coil(Magnetic magnetic, OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold = defaults.harmonicAmplitudeThreshold, std::optional<std::vector<int8_t>> customCurrentDirectionPerWinding = std::nullopt, std::optional<CoilMesherModels> coilMesherModel = std::nullopt);
    std::vector<Field> generate_mesh_induced_coil(Magnetic magnetic, OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold = defaults.harmonicAmplitudeThreshold);
    std::vector<size_t> get_common_harmonic_indexes(OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold);
    // meshAllWindows: mesh every distinct winding-window region (multi-column cores)
    // instead of only window 0. Painting wants it; the leakage energy integrator must
    // NOT use it — its revolution bookkeeping already accounts for the full turn from
    // one window's cross-section, so meshing both sides would double-count energy.
    static std::pair<Field, double> generate_mesh_induced_grid(Magnetic magnetic, double frequency, size_t numberPointsX, size_t numberPointsY, bool ignoreTurns = false, bool includeInsideTurns = true, bool meshAllWindows = false);
    // Uniform cell-centred grid over the WHOLE core winding window of a single-window core with a rectangular
    // window: x from the main column face to the lateral column face, y over the window height, conductors
    // included. The leakage energy integral needs every point where the field stores energy: the copper, the
    // space between separated sections, and the bobbin walls and clearances beyond the bobbin window
    // (ABT #1240). Throws for toroids and multi-window cores, whose windows this frame does not describe.
    static std::pair<Field, double> generate_mesh_core_winding_window_grid(Magnetic magnetic, double frequency, size_t numberPointsX, size_t numberPointsY);
};

class CoilMesherModel {
  private:
  public:
    std::string method_name = "Default";
    // Reference parameters: these run once PER TURN inside the adviser hot path; the old
    // by-value Turn/Wire/Core signatures deep-copied the whole Core (including the cached
    // material datasets) for every turn. Wire/Core are NON-const references only because
    // their getters (get_maximum_conducting_*, get_initial_permeability) sit on non-const
    // resolve chains — the models never mutate them.
    virtual std::vector<FieldPoint> generate_mesh_inducing_turn(const Turn& turn, Wire& wire, std::optional<size_t> turnIndex, std::optional<double> turnLength, Core& core) = 0;
    // core: optional (nullptr allowed), needed only to tell a genuinely lateral-column
    // crossing (rectangular multi-column winding, ABT #227.2) apart from a toroidal
    // turn's own outer-return crossing, which already gets dedicated Kelvin-image
    // handling on the inducing side and must not also pick up a second induced-side
    // sample here.
    virtual std::vector<FieldPoint> generate_mesh_induced_turn(const Turn& turn, Wire& wire, std::optional<size_t> turnIndex = std::nullopt, const Core* core = nullptr) = 0;
    static std::shared_ptr<CoilMesherModel> factory(CoilMesherModels modelName);

};

class CoilMesherCenterModel : public CoilMesherModel {
  public:
    std::vector<FieldPoint> generate_mesh_inducing_turn(const Turn& turn, [[maybe_unused]] Wire& wire, std::optional<size_t> turnIndex, std::optional<double> turnLength, Core& core);
    std::vector<FieldPoint> generate_mesh_induced_turn(const Turn& turn, [[maybe_unused]] Wire& wire, std::optional<size_t> turnIndex = std::nullopt, const Core* core = nullptr);
};

// // Based on Improved Analytical Calculation of High Frequency Winding Losses in Planar Inductors by Xiaohui Wang
// // https://sci-hub.wf/10.1109/ECCE.2018.8558397
class CoilMesherWangModel : public CoilMesherModel {
  public:
    std::vector<FieldPoint> generate_mesh_induced_turn(const Turn& turn, Wire& wire, std::optional<size_t> turnIndex = std::nullopt, const Core* core = nullptr);
    std::vector<FieldPoint> generate_mesh_inducing_turn(const Turn& turn, Wire& wire, std::optional<size_t> turnIndex, std::optional<double> turnLength, [[maybe_unused]] Core& core);
};



} // namespace OpenMagnetics