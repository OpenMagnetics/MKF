#pragma once
#include "MAS.hpp"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Wire.h"
#include "support/Utils.h"

#include <cmath>
#include <complex>
#include <optional>
#include <string>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <streambuf>
#include <string_view>
#include <vector>

using namespace MAS;

namespace OpenMagnetics {

// Labels of the toroid inducing points (CoilMesherCenterModel::generate_mesh_inducing_turn). The ring splits the plane
// into the bore and the outside; each point acts only in its region (MagneticField). A "reply" is the core's answer to
// a conductor (its Kelvin image, or the centre line current that carries a bore crossing's MMF outside the core).
inline constexpr std::string_view toroidBoreConductorLabel = "toroid bore conductor";
inline constexpr std::string_view toroidBoreReplyLabel = "toroid bore reply";
inline constexpr std::string_view toroidExteriorConductorLabel = "toroid exterior conductor";
inline constexpr std::string_view toroidExteriorReplyLabel = "toroid exterior reply";


enum class CoilMesherModels : int {
    WANG,
    CENTER
};

// ABT #1304: a uniform current sheet on a straight segment of the winding-window boundary, per unit depth. It is
// the equivalent source of an MMF drop on that boundary -- the mouth of a gap, or the ferrite along a wall -- and
// every lattice image of it is returned as its own sheet with its image factor already applied.
struct WindowWallSheet {
    std::vector<double> start;
    std::vector<double> end;
    // Total current of the sheet per ampere-turn of NET winding MMF in the window, sign included. The sheets of one
    // window (images excluded) sum to exactly -1: they return the winding's MMF, which is what makes each lattice
    // cell neutral.
    // Multi-window cores (ABT #2267): per ampere-turn of the MMF of winding `winding` (N_w I_w, its polarity included),
    // and the sheets of window `window` sum to minus the current that winding encloses in that window.
    double current;
    std::optional<size_t> winding = std::nullopt;
    size_t window = 0;
};

// Labels of the CENTER mesher's points: "image" marks the core's reply to a conductor (an image), and on a multi-window
// core (ABT #2267) "@k" the window k the point belongs to ("image@k" an image in it). A point acts only in its window.
inline std::string window_label(bool image, std::optional<size_t> window) {
    return std::string(image ? "image" : "") + (window ? "@" + std::to_string(window.value()) : std::string());
}
inline bool is_image_label(const std::optional<std::string>& label) {
    return label && label->rfind("image", 0) == 0;
}
inline std::optional<size_t> window_of_label(const std::optional<std::string>& label) {
    if (!label) {
        return std::nullopt;
    }
    const auto at = label->find('@');
    if (at == std::string::npos) {
        return std::nullopt;
    }
    return static_cast<size_t>(std::stoul(label->substr(at + 1)));
}

// A winding-window region of a multi-window core (ABT #2267): the air between two x-adjacent legs, the frame its
// conductors and wall sheets are imaged in.
struct WindowRegion {
    double x0;
    double y0;
    double x1;
    double y1;
    bool contains(double x, double y) const { return x > x0 && x < x1 && y > y0 && y < y1; }
};


// The inducing mesh of a coil together with the PHASE of every winding's current at every
// meshed harmonic (MAS excitation convention, 2026-09-24).
//
// fieldPerHarmonic carries, per turn, the signed peak AMPLITUDE of the harmonic current
// (MKF's MAS harmonic amplitude x the turn's current divider x the winding's direction c_k),
// exactly as generate_mesh_inducing_coil always did. currentPhasePerHarmonicPerWinding[i][k]
// is the phase (rad) of winding k's current at the harmonic of fieldPerHarmonic[i], taken
// from the winding's Fourier coefficient at that harmonic (the one its MAS amplitude comes
// from, ABT #1880) and referred to the gauge winding (the
// reference winding, see calculate_current_phase_per_winding), so the complex current of a
// turn of winding k is  value * exp(j * phase[k]).  Inducing points that belong to no turn
// (equivalent fringing sources) carry the gauge phase 0.
struct InducingCoilMesh {
    std::vector<Field> fieldPerHarmonic;
    std::vector<std::vector<double>> currentPhasePerHarmonicPerWinding;
    // ABSOLUTE phase (rad, on the waveforms' own time base) of the gauge at the harmonic of
    // fieldPerHarmonic[i]: subtracting it refers any other signal on the same time base (the
    // magnetizing current, whose phase drives the gap fringing field) to the same gauge.
    // Empty when no winding carries that harmonic.
    std::vector<std::optional<double>> gaugePhasePerHarmonic;
    // ABT #1314: the peak current per winding each harmonic's points were scaled by (before divider and direction). With
    // unitCurrentForSilentParallels, a winding of several parallels that carries nothing at a harmonic is meshed at 1 A so
    // a parallel current split can place its circulating current there; every consumer then rescales by this amplitude.
    std::vector<std::vector<double>> amplitudePerHarmonicPerWinding;
    // The peak current each winding carries at each harmonic (calculate_current_phase_per_winding): its MAS amplitude, 0
    // where it carries none (ABT #1880: a harmonic zero both in the list and in the waveform). No silent-parallel 1 A.
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
    // harmonic indexes: arg of the waveform's Fourier coefficient at the frequency MAS lists for
    // that harmonic, the same coefficient its MAS amplitude comes from (ABT #1880). Referred to a gauge: the reference winding's phase at that
    // harmonic, or — when the reference winding carries none of that harmonic — the first
    // winding (by index) that does. A global rotation leaves every |H|^2 unchanged; the gauge
    // only fixes WHICH component is reported as in-phase, and makes a single winding and
    // windings in exact antiphase (all phases 0) reproduce the amplitude-only field exactly.
    // A winding that carries no current at a harmonic has no phase there: its MAS amplitude is
    // zero, or the harmonic is zero both in the list and in its waveform
    // (WaveformProcessor::is_zero_harmonic). Its entry is a placeholder 0 that must not be used:
    // its current is 0 (currentAmplitudePerHarmonic), so it contributes nothing.
    // Throws when an excited winding has no waveform or its waveform does not contain the
    // harmonic MAS lists.
    // Result: [position in harmonicIndexes][winding index].
    // gaugePhasePerHarmonic, when given, receives the gauge's absolute phase per harmonic.
    // currentAmplitudePerHarmonic, when given, receives [position][winding] the peak current each
    // winding carries there: its MAS amplitude, 0 where it carries none.
    static std::vector<std::vector<double>> calculate_current_phase_per_winding(Coil coil, OperatingPoint operatingPoint, const std::vector<size_t>& harmonicIndexes, std::vector<std::optional<double>>* gaugePhasePerHarmonic = nullptr, std::vector<std::vector<double>>* currentAmplitudePerHarmonic = nullptr);
    // Peak phasor X of a waveform at harmonicFrequency, x(t) = Re(X exp(+j w t)): 2 c_k from the
    // same Fourier coefficients its MAS amplitude comes from (WaveformProcessor::calculate_harmonic_phasor,
    // ABT #1880), so |X| is that amplitude and the phase is consistent with it by construction.
    static std::complex<double> calculate_harmonic_phasor(const Waveform& waveform, double excitationFrequency, double harmonicFrequency);
    InducingCoilMesh generate_mesh_inducing_coil_phasors(Magnetic magnetic, OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold = defaults.harmonicAmplitudeThreshold, std::optional<std::vector<int8_t>> customCurrentDirectionPerWinding = std::nullopt, std::optional<CoilMesherModels> coilMesherModel = std::nullopt, bool unitCurrentForSilentParallels = false);
    // Amplitude-and-direction part only (generate_mesh_inducing_coil_phasors().fieldPerHarmonic).
    // The phases are NOT in these fields: a caller that sums them treats every winding as in
    // phase with the reference. Field computations must use the phasor variant.
    std::vector<Field> generate_mesh_inducing_coil(Magnetic magnetic, OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold = defaults.harmonicAmplitudeThreshold, std::optional<std::vector<int8_t>> customCurrentDirectionPerWinding = std::nullopt, std::optional<CoilMesherModels> coilMesherModel = std::nullopt);
    std::vector<Field> generate_mesh_induced_coil(Magnetic magnetic, OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold = defaults.harmonicAmplitudeThreshold);
    std::vector<size_t> get_common_harmonic_indexes(OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold);
    // ABT #1304: the magnetic circuit's MMF drops as sheets on the boundary of a single rectangular winding window,
    // with their complete-cell image lattice. Each gap's share is its reluctance over the total (MKF's own
    // composition: main-column gaps in series, lateral branches in parallel); the ferrite's share is spread along
    // the four walls in proportion to length x tangential H, H = B / mu with B from the flux through that wall's
    // limb. Throws for toroids and multi-window cores, whose boundary this frame does not describe.
    static std::vector<WindowWallSheet> generate_window_wall_sheets(Magnetic magnetic, std::optional<OperatingPoint> operatingPoint);
    // The same sources before imaging: one sheet per gap mouth and per wall's ferrite share, currents per ampere-turn
    // of net MMF summing to -1 (AxisymmetricWindow uses them as the window's boundary data).
    static std::vector<WindowWallSheet> generate_window_wall_sheet_sources(Magnetic magnetic, std::optional<OperatingPoint> operatingPoint);
    // Multi-window cores (ABT #2267): the window regions between x-adjacent legs, in leg order, and the wall sheet
    // sources per winding and window from the column reluctance network (ReluctanceNetwork): each gap's drop is its
    // reluctance times its column's flux, each leg face carries the leg's ferrite drop and the two plates the yoke
    // segment's, so the sources of window k for winding w sum to minus the current w encloses there.
    static std::vector<WindowRegion> window_regions(Core core);
    // A gap's slot mouth on the window wall carries the gap's MMF drop as slotMouthPieces sheets of equal current, at the
    // boundaries this returns: fractions of the mouth from its leg corner (0) to its far end (1). A mouth inside the wall
    // (atWindowCorner false) is the slot between two leg corners, opening into the half-space; one at the window's end
    // (atWindowCorner true) is the slot between a leg corner and a plate face flush with the window's end wall.
    static constexpr size_t slotMouthPieces = 8;
    static const std::vector<double>& slot_mouth_boundaries(bool atWindowCorner);

    // ABT #2384: a gap's slot as a conductor near its mouth sees it. The lattice images a conductor in the column face as
    // if the face were ferrite across the mouth too; the slot lets the conductor's own field into it. The map of the slot
    // (gap length g) onto the upper half t-plane with the irons on the real axis and the leg corner at t = 1,
    //     z = A [s + j ln((1 + j s) / t)],  s = sqrt(t - 1) sqrt(t + 1)   (two leg corners, dz/dt = A s / t),
    //     z = A [2 r - 2 atan r],           r = sqrt(t - 1)               (a leg corner and the window's end wall),
    // A = -j g / pi, z measured from the leg corner with the window at Re z > 0 and the slot along 0 < Im z < g, carries
    // the exact field of a line current I at t0: the same current at conj(t0) is its image. That field differs from the
    // lattice's by the field of a sheet on the mouth (with the coincident image the lattice gives it) whose current
    // between two of its points is the change of I (arg(t - t0) + arg(t - conj t0)) / (2 pi) between them: within 0.2-0.7 %
    // of the exact field on 32 pieces (scratch check against the map's own field, 2026-11-15). The legs stay at one
    // potential, so a mouth's pieces carry no net current.
    struct SlotMouth {
        double wallX = 0;
        // The leg corner the map puts at t = 1, and the gap's length along the wall.
        double cornerY = 0;
        double width = 0;
        // +1 when the window is at x > wallX (a column face on the window's left), -1 when it is at x < wallX.
        double towardsWindow = 1;
        // +1 when the slot runs from cornerY to cornerY + width, -1 when it runs to cornerY - width.
        double alongSlot = 1;
        // The far side of the mouth is a plate face flush with the window's end wall, not a second leg corner.
        bool atWindowCorner = false;
        // The pieces' ends from the corner (0) to the far side (width), and their images in the t-plane.
        std::vector<double> boundaryOffsets;
        std::vector<std::complex<double>> boundaryT;
    };
    static constexpr size_t slotCorrectionPieces = 32;
    // The mouths of a single-window core's gaps on its two walls, as generate_window_wall_sheet_sources opens them.
    static std::vector<SlotMouth> slot_mouths(Magnetic magnetic);
    // The map of a slot of unit width, and its inverse (Newton, kept in the upper half plane; throws if it does not converge).
    // `start`, the preimage of a nearby point, is tried before the far-field starts.
    static std::complex<double> slot_map(std::complex<double> t, bool atWindowCorner);
    static std::complex<double> slot_map_inverse(std::complex<double> z, bool atWindowCorner, std::optional<std::complex<double>> start = std::nullopt);
    // Per piece of the mouth, its sheet current per ampere of a line current at (x, y) in the window. `preimage` (in/out):
    // a nearby point's preimage to start from, replaced by this point's.
    static std::vector<double> slot_correction_currents(const SlotMouth& mouth, double x, double y, std::optional<std::complex<double>>* preimage = nullptr);
    // The mouth's pieces as sheets of unit current, before imaging, in piece order.
    static std::vector<WindowWallSheet> slot_correction_pieces(const SlotMouth& mouth);
    // Sheets of a single-window core imaged in its complete-cell lattice, as generate_window_wall_sheets images the walls.
    static std::vector<WindowWallSheet> image_window_wall_sheets(const std::vector<WindowWallSheet>& originals, Core core);
    // Appends the mouth of a gap on the wall x = wallX, from yStart to yEnd, as those pieces. cornerAtStart says which end
    // is the leg corner when the mouth is at a window corner.
    static void append_slot_mouth(std::vector<WindowWallSheet>& sheets, double wallX, double yStart, double yEnd, double current, bool atWindowCorner,
                                  bool cornerAtStart, std::optional<size_t> winding = std::nullopt, size_t window = 0);
    static std::vector<WindowWallSheet> generate_multi_window_wall_sheet_sources(Magnetic magnetic, std::optional<OperatingPoint> operatingPoint);
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
    // A round-post window seen as a solid of revolution (AxisymmetricWindow::applies, decided once per coil by
    // CoilMesher::generate_mesh_inducing_coil on the neutral path): turns imaged in the plates only.
    std::optional<bool> axisymmetricWindow;
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
