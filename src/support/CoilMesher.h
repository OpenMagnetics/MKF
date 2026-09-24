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


// The inducing mesh of a coil together with the PHASE of every winding's current at every
// meshed harmonic (MAS excitation convention, 2026-09-24).
//
// fieldPerHarmonic carries, per turn, the signed peak AMPLITUDE of the harmonic current
// (MKF's MAS harmonic amplitude x the turn's current divider x the winding's direction c_k),
// exactly as generate_mesh_inducing_coil always did. currentPhasePerHarmonicPerWinding[i][k]
// is the phase (rad) of winding k's current at the harmonic of fieldPerHarmonic[i], taken
// from a DFT of the winding's current waveform and referred to the gauge winding (the
// reference winding, see calculate_current_phase_per_winding), so the complex current of a
// turn of winding k is  value * exp(j * phase[k]).  Inducing points that belong to no turn
// (equivalent fringing sources) carry the gauge phase 0.
struct InducingCoilMesh {
    std::vector<Field> fieldPerHarmonic;
    std::vector<std::vector<double>> currentPhasePerHarmonicPerWinding;
    // MAS harmonic amplitude (peak, A) of each winding's current at the harmonic of
    // fieldPerHarmonic[i]; 0 when the winding lists none. With the phases and the directions it
    // gives the magnetizing current phasor i_m = sum_k c_k N_k I_k / N_r (the gap fringing source).
    std::vector<std::vector<double>> currentAmplitudePerHarmonicPerWinding;
};

class CoilMesher {
  private:
  protected:
    double _quickModeForManyHarmonicsThreshold = 1;
  public:
    // MAS excitation convention (2026-09-24). The reference winding r is the FIRST winding
    // whose isolationSide is "primary"; throws if the coil has none.
    static size_t get_reference_winding_index(Coil coil);
    // MAS excitation convention: primary-side windings are PASSIVE (+ into the dot), every
    // other winding SOURCES (+ out of the dot), so the physical dot-referenced current of
    // winding k is c_k * i_k with c_k = +1 for isolationSide "primary" and -1 otherwise.
    static std::vector<int8_t> calculate_current_direction_per_winding(Coil coil);
    // Phase (rad, e^{+j w t} convention) of each winding's current at each of the given
    // harmonic indexes, from a DFT of the winding's current waveform at the frequency MAS
    // lists for that harmonic. Referred to a gauge: the reference winding's phase at that
    // harmonic, or — when the reference winding carries none of that harmonic — the first
    // winding (by index) that does. A global rotation leaves every |H|^2 unchanged; the gauge
    // only fixes WHICH component is reported as in-phase, and makes a single winding and
    // windings in exact antiphase (all phases 0) reproduce the amplitude-only field exactly.
    // A winding whose MAS amplitude at a harmonic is zero gets phase 0 (it contributes nothing).
    // Throws when an excited winding has no waveform or its waveform does not contain the
    // harmonic MAS lists.
    // Result: [position in harmonicIndexes][winding index].
    static std::vector<std::vector<double>> calculate_current_phase_per_winding(Coil coil, OperatingPoint operatingPoint, const std::vector<size_t>& harmonicIndexes);
    InducingCoilMesh generate_mesh_inducing_coil_phasors(Magnetic magnetic, OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold = defaults.harmonicAmplitudeThreshold, std::optional<std::vector<int8_t>> customCurrentDirectionPerWinding = std::nullopt, std::optional<CoilMesherModels> coilMesherModel = std::nullopt);
    // Amplitude-and-direction part only (generate_mesh_inducing_coil_phasors().fieldPerHarmonic).
    // The phases are NOT in these fields: a caller that sums them treats every winding as in
    // phase with the reference. Field computations must use the phasor variant.
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