#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "constructive_models/Coil.h"
#include "constructive_models/Core.h"
#include "support/CoilMesher.h"

namespace OpenMagnetics {

// The winding window of a round-post core seen as a solid of revolution: a ring-shaped box a < r < b,
// zBottom < z < zTop, bounded by ferrite (mu -> inf) on all four sides. Used by the neutral-window field
// (MMF_BOUNDARY on BINNS_LAWRENSON) on the cores OMFEM also solves axisymmetric (see applies()).
//
// A turn is a ring, not a line: its field is the exact ring field (complete elliptic integrals). The plates
// (z = zBottom, zTop) are planes, so a ring's image in them is again a ring, exact. The post face and the outer
// wall are cylinders, which have no ring image; their condition (tangential H = 0, or the MMF a gap or the ferrite
// drops there) is met by a correction field regular in the window: a cosine series in z (zero H_r on the plates)
// with the modified Bessel functions I0, K0 (H_z) and I1, K1 (H_r) in r. A ring's own H_z on either cylinder,
// together with all its plate images, has closed-form cosine coefficients, so the correction needs no sampling:
//     on r = a < r':  H_z,n = (2 I r' / B) k I0(k a) K1(k r') cos(k (z' - zBottom)),  H_z,0 = I / B
//     on r = b > r':  H_z,n = -(2 I r' / B) k K0(k b) I1(k r') cos(k (z' - zBottom)), H_z,0 = 0
// with k = n pi / B and B = zTop - zBottom. Checked against OMFEM's axisymmetric solve of an ETD 29 flyback
// (primary alone, 2 kHz): within a few per cent across the window, where the planar section read 1.13-1.47x.
class AxisymmetricWindow {
public:
    // Round main column, one winding window, not a toroid, and a functional (subtractive or additive) gap or at most
    // one lateral limb: the cores whose window OMFEM's 2D solve treats as a solid of revolution. An ungapped round
    // post between two discrete outer legs stays planar, as OMFEM solves it.
    static bool applies(const Core& core);
    // The same for a wound magnetic: the core's rule, and every turn's section inside the window. A stated layout that
    // overflows it (laid out as the design gives it, ABT #2096) has copper in the ferrite, where the correction is not
    // defined, and stays on the planar section.
    static bool applies(const Core& core, const Coil& coil);
    static AxisymmetricWindow from_core(const Core& core);

    AxisymmetricWindow(double postRadius, double outerRadius, double zBottom, double zTop);

    double post_radius() const { return _a; }
    double outer_radius() const { return _b; }
    double bottom() const { return _zBottom; }
    double top() const { return _zTop; }
    bool contains(double r, double z) const;

    // (H_r, H_z) of a ring of current I (in the sense of MKF's planar line current, + phi) at (rs, zs), at (r, z).
    static std::pair<double, double> ring_field(double current, double rs, double zs, double r, double z);
    // The same ring with its images in the two plates (planeImages reflections each way).
    std::pair<double, double> ring_field_with_plate_images(double current, double rs, double zs, double r, double z) const;
    // H of a ring minus that of MKF's planar line current through the same section point: the curvature part,
    // added to a planar rectangle's exact field for a flat conductor. The distance is held at minimumDistance or more.
    static std::pair<double, double> ring_minus_line_field(double current, double rs, double zs, double r, double z, double minimumDistance);

    // H_z wanted on the two cylinders minus what the sources give there, as cosine coefficients in z.
    struct WallMismatch {
        std::vector<double> post;
        std::vector<double> outer;
        double postMean = 0;
        double outerMean = 0;
    };
    WallMismatch empty_mismatch() const;
    // A ring source and its plate images: its H_z on the cylinders is to be cancelled.
    void add_ring(WallMismatch& mismatch, double current, double rs, double zs) const;
    // A tangential H wanted on a band of a cylinder (a gap mouth, or the ferrite's share of the MMF).
    void add_wall_field(WallMismatch& mismatch, bool postFace, double zStart, double zEnd, double tangentialField) const;

    struct Correction {
        std::vector<double> p;
        std::vector<double> q;
        double uniform = 0;
        size_t modes = 0;  // the last mode that carries anything (solve sets it)
    };
    Correction solve(const WallMismatch& mismatch) const;
    // (H_r, H_z) of the correction at a point of the window.
    std::pair<double, double> correction_field(const Correction& correction, double r, double z) const;

    // The correction's basis at one point, for evaluating several corrections there.
    struct PointBasis {
        size_t modes = 0;
        std::vector<double> radialP;
        std::vector<double> radialQ;
        std::vector<double> axialP;
        std::vector<double> axialQ;
    };
    PointBasis basis_at(double r, double z) const;
    static std::pair<double, double> correction_field(const Correction& correction, const PointBasis& basis);
    // The same correction split by the wall it answers: the post face's reply (the K series, decaying away from r = a)
    // and the outer wall's (the I series, decaying away from r = b), the uniform mode shared evenly. A foil between a
    // point and one wall screens that wall's part only (what the planar section's image in that wall would cross).
    struct WallParts {
        std::pair<double, double> post;
        std::pair<double, double> outer;
    };
    static WallParts correction_field_by_wall(const Correction& correction, const PointBasis& basis);

    // The neutral window's wall sources (CoilMesher::generate_window_wall_sheet_sources, per ampere-turn of net MMF) in
    // this window: the gaps' and the ferrite's MMF on the two cylinders is the tangential H the correction meets there
    // (just inside a wall MKF's planar sheets give -I/L on the post face and +I/L on the outer wall, I a sheet's current);
    // the plates' share stays a source, rings on the plates (Gauss-Legendre, 8 per sheet) doubled by their plate images.
    struct PlateRing {
        double radius;
        double height;
        double current;
    };
    // A band of tangential H the correction meets on a cylinder (a gap mouth piece, or the ferrite's share).
    struct WallBand {
        bool postFace;
        double zStart;
        double zEnd;
        double tangentialField;
    };
    struct WallSources {
        std::vector<PlateRing> plateRings;
        Correction correction;
        std::vector<WallBand> bands;
    };
    WallSources wall_sources(const std::vector<WindowWallSheet>& sheetSources) const;
    // The bands' modes past the ones the correction keeps, at a point nearer a cylinder than they resolve: (H_r, H_z) of
    // the post's and of the outer wall's, in the planar limit the high modes reach (see the definition).
    struct BandTail {
        std::pair<double, double> post{0.0, 0.0};
        std::pair<double, double> outer{0.0, 0.0};
    };
    BandTail band_tail(const WallSources& sources, double r, double z) const;
    // (Hx, Hy) of the wall sources per ampere-turn of net MMF at a section point (x the radius; x < 0 is the far side
    // of the axis). The correction acts inside the window only.
    std::pair<double, double> wall_field(const WallSources& sources, double pointX, double pointY) const;
    // The same, one part per source and per wall, for foil screening: each plate ring (at its own position) and the
    // correction's post and outer parts (at the wall, the point's height), as section fields (x < 0 mirrored).
    struct WallFieldPart {
        double sourceX;
        double sourceY;
        double fieldX;
        double fieldY;
    };
    std::vector<WallFieldPart> wall_field_parts(const WallSources& sources, double pointX, double pointY) const;

    // Flux through the circle of radius r at height z of a unit-current ring at (rs, zs): 2 pi r A_phi, in Wb/A.
    static double ring_flux(double rs, double zs, double r, double z);
    double ring_flux_with_plate_images(double rs, double zs, double r, double z) const;
    // Flux through the circle (r, z) of a correction field, up to a constant (the same for every circle): 2 pi mu0 psi.
    double correction_flux(const Correction& correction, double r, double z) const;
    double wall_flux(const WallSources& sources, double r, double z) const;
    // correction_flux at (r, z) minus that at a reference circle, per unit coefficient (every mode kept): one basis per
    // target, dotted with each source's correction, where the split reads many corrections at one target.
    struct FluxBasis {
        double uniform = 0;
        std::vector<double> p;
        std::vector<double> q;
    };
    FluxBasis flux_basis_between(double r, double z, double referenceR, double referenceZ) const;
    static double correction_flux(const Correction& correction, const FluxBasis& basis);

    size_t modes() const { return _modes; }
    static constexpr int planeImages = 4;

private:
    double _a;
    double _b;
    double _zBottom;
    double _zTop;
    double _height;
    size_t _modes;
    std::vector<double> _k;
    // Per mode: I0(ka)/I0(kb) and K0(kb)/K0(ka), and the 2x2 determinant they form.
    std::vector<double> _postRatio;
    std::vector<double> _outerRatio;
    // Per mode: 1 / (I0(kb) e^-kb) and 1 / (K0(ka) e^ka), the basis's normalisations.
    std::vector<double> _inverseScaledI0Outer;
    std::vector<double> _inverseScaledK0Post;
    // Per mode: I0(ka) e^-ka and K0(kb) e^kb, a ring's wall coefficients' fixed factors.
    std::vector<double> _scaledI0Post;
    std::vector<double> _scaledK0Outer;
};

}  // namespace OpenMagnetics
