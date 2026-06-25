#include "physical_models/WindingLossesPEEC.h"
#include "physical_models/Resistivity.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "constructive_models/Wire.h"
#include "constructive_models/Core.h"
#include "support/CoilMesher.h"
#include "support/Settings.h"
#include "support/Exceptions.h"
#include "Definitions.h"

#include <Eigen/Dense>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <numbers>

namespace OpenMagnetics {

using cd = std::complex<double>;
static const double MU0 = 4 * std::numbers::pi * 1e-7;

// Geometric mean distance of a rectangle to itself (Rosa), used for the
// self partial inductance term so ln(distance) doesn't blow up at d=0.
static double self_gmd(double w, double h) {
    return 0.2235 * (w + h);
}

// Average of ln(distance) over two rectangular cells = the 2D mutual
// partial-inductance kernel <ln d>_ij. Uses 3-point Gauss-Legendre per
// dimension (81 pairs). For distinct (non-overlapping) cells the integrand is
// smooth at the interior Gauss points, so this converges as the mesh refines —
// which the centroid point-log approximation does NOT (that was the source of
// the mesh non-convergence).
static double avg_ln_between_cells(double xi, double yi, double wi, double hi,
                                   double xj, double yj, double wj, double hj) {
    static const double gp[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
    static const double gw[3] = {0.5555555555555556, 0.8888888888888889, 0.5555555555555556};
    double acc = 0, wsum = 0;
    for (int ax = 0; ax < 3; ++ax) {
        for (int ay = 0; ay < 3; ++ay) {
            double pax = xi + 0.5 * wi * gp[ax];
            double pay = yi + 0.5 * hi * gp[ay];
            double wA = gw[ax] * gw[ay];
            for (int bx = 0; bx < 3; ++bx) {
                for (int by = 0; by < 3; ++by) {
                    double pbx = xj + 0.5 * wj * gp[bx];
                    double pby = yj + 0.5 * hj * gp[by];
                    double wB = gw[bx] * gw[by];
                    double dx = pax - pbx, dy = pay - pby;
                    double d = std::sqrt(dx * dx + dy * dy);
                    acc += wA * wB * std::log(std::max(d, 1e-15));
                    wsum += wA * wB;
                }
            }
        }
    }
    return acc / wsum;
}

// Symmetric GRADED cell edges on [-D/2, D/2]: smallest cell ~ cellMin at both
// surfaces (where current crowds at high frequency), growing geometrically by
// `ratio` toward the centre. Resolves the skin depth with far fewer filaments
// than a uniform grid -> small, well-conditioned matrices. Returns edges
// (size = nCells+1), normalised to span exactly D.
static std::vector<double> graded_edges(double D, double cellMin, double ratio, size_t maxCellsPerSide) {
    cellMin = std::min(cellMin, D / 2);
    if (cellMin <= 0) cellMin = D / 2;
    std::vector<double> half; // surface -> centre, increasing widths
    double c = cellMin, sum = 0;
    while (sum + c < D / 2 && half.size() + 1 < maxCellsPerSide) {
        half.push_back(c);
        sum += c;
        c *= ratio;
    }
    half.push_back(std::max(D / 2 - sum, cellMin * 0.25)); // remaining centre half-cell
    std::vector<double> widths(half); // surface -> centre ...
    for (auto it = half.rbegin(); it != half.rend(); ++it) widths.push_back(*it); // ... centre -> surface
    double total = 0;
    for (double w : widths) total += w;
    double scale = D / total;
    std::vector<double> edges;
    edges.reserve(widths.size() + 1);
    double x = -D / 2;
    edges.push_back(x);
    for (double w : widths) {
        x += w * scale;
        edges.push_back(x);
    }
    return edges;
}

// Subdivide a turn's cross-section into a graded tensor mesh, refined toward the
// surfaces in BOTH dimensions. cellMin ~ skinDepth/3 (penetration resolution);
// maxCellsPerSide bounds the per-dimension cell count.
static void mesh_turn_into_filaments(const Turn& turn, Wire wire, size_t turnIndex, size_t windingIndex,
                                     double cellMin, double ratio, size_t maxCellsPerSide,
                                     std::vector<WindingLossesPEEC::Filament>& out) {
    double cx = turn.get_coordinates()[0];
    double cy = turn.get_coordinates()[1];
    double turnLength = turn.get_length();

    bool isRound = (wire.get_type() == WireType::ROUND);
    double fullW, fullH, radius = 0;
    if (isRound) {
        radius = resolve_dimensional_values(wire.get_conducting_diameter().value()) / 2;
        fullW = fullH = 2 * radius;
    }
    else {
        fullW = wire.get_maximum_conducting_width();
        fullH = wire.get_maximum_conducting_height();
    }

    auto ex = graded_edges(fullW, cellMin, ratio, maxCellsPerSide);
    auto ey = graded_edges(fullH, cellMin, ratio, maxCellsPerSide);

    struct Cell { double x, y, w, h; };
    std::vector<Cell> cells;
    double boxAreaInside = 0;
    for (size_t i = 0; i + 1 < ex.size(); ++i) {
        for (size_t j = 0; j + 1 < ey.size(); ++j) {
            double w = ex[i + 1] - ex[i];
            double h = ey[j + 1] - ey[j];
            double fx = cx + 0.5 * (ex[i] + ex[i + 1]);
            double fy = cy + 0.5 * (ey[j] + ey[j + 1]);
            if (isRound) {
                double rx = fx - cx, ry = fy - cy;
                if (rx * rx + ry * ry > radius * radius) continue;
            }
            cells.push_back({fx, fy, w, h});
            boxAreaInside += w * h;
        }
    }
    if (cells.empty()) return;

    double conductingArea = isRound ? (std::numbers::pi * radius * radius) : (fullW * fullH);
    // conserve total conducting area (exact for rectangles; rescales the round
    // bounding-box tiling back to the true circle area).
    double areaScale = (boxAreaInside > 0) ? conductingArea / boxAreaInside : 1.0;

    for (auto& cl : cells) {
        WindingLossesPEEC::Filament f;
        f.x = cl.x;
        f.y = cl.y;
        f.area = cl.w * cl.h * areaScale;
        f.width = cl.w;
        f.height = cl.h;
        f.turnIndex = turnIndex;
        f.windingIndex = windingIndex;
        f.turnLength = turnLength;
        out.push_back(f);
    }
}

WindingLossesOutput WindingLossesPEEC::calculate_losses(Magnetic magnetic, OperatingPoint operatingPoint, double temperature) {
    auto& settings = Settings::GetInstance();
    (void)settings;

    // --- core geometry (self-heal unprocessed cores, like MagneticField) ---
    auto core = magnetic.get_core();
    if (!core.get_processed_description()) {
        core.process_data();
    }
    if (!core.is_gap_processed()) {
        core.process_gap();
    }
    magnetic.set_core(core);

    if (core.get_shape_family() == CoreShapeFamily::T) {
        throw NotImplementedException("WindingLossesPEEC: toroidal cores not supported yet");
    }

    auto coil = magnetic.get_coil();
    if (!coil.get_turns_description()) {
        throw CoilNotProcessedException("WindingLossesPEEC: coil has no turns description");
    }
    auto turns = coil.get_turns_description().value();
    size_t numberTurns = turns.size();

    // --- highest significant frequency (drives the mesh resolution) ---
    auto excitationsForMesh = operatingPoint.get_excitations_per_winding();
    double maxSignificantFrequency = 0;
    double globalMaxAmp = 0;
    for (auto& exc : excitationsForMesh) {
        if (exc.get_current() && exc.get_current()->get_harmonics()) {
            for (double a : exc.get_current()->get_harmonics()->get_amplitudes()) {
                globalMaxAmp = std::max(globalMaxAmp, std::fabs(a));
            }
        }
    }
    for (auto& exc : excitationsForMesh) {
        if (exc.get_current() && exc.get_current()->get_harmonics()) {
            auto amps = exc.get_current()->get_harmonics()->get_amplitudes();
            auto freqs = exc.get_current()->get_harmonics()->get_frequencies();
            for (size_t h = 0; h < amps.size() && h < freqs.size(); ++h) {
                if (globalMaxAmp > 0 && std::fabs(amps[h]) / globalMaxAmp >= _harmonicAmplitudeThreshold) {
                    maxSignificantFrequency = std::max(maxSignificantFrequency, freqs[h]);
                }
            }
        }
    }
    if (maxSignificantFrequency <= 0) {
        maxSignificantFrequency = 1.0;
    }

    // --- mesh every turn into filaments; reject litz (impractical to mesh strands) ---
    // Resolution is chosen per turn so the filament pitch ~ skinDepth/3 in each
    // direction at the highest significant frequency.
    std::vector<Filament> filaments;
    for (size_t turnIndex = 0; turnIndex < numberTurns; ++turnIndex) {
        auto turn = turns[turnIndex];
        size_t windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        auto wire = coil.resolve_wire(windingIndex);
        if (wire.get_type() == WireType::LITZ) {
            throw NotImplementedException("WindingLossesPEEC: litz must use the analytical model (strand meshing is impractical); homogenise the bundle to enable PEEC");
        }
        double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, maxSignificantFrequency, temperature);
        // graded mesh: ~skinDepth/factor cells at the surfaces, growing by ratio.
        double cellMin = skinDepth / _cellMinFactor;
        size_t maxCellsPerSide = (numberTurns > 30 && _maxCellsPerSide > 8) ? 8 : _maxCellsPerSide;
        mesh_turn_into_filaments(turn, wire, turnIndex, windingIndex, cellMin, _gradingRatio, maxCellsPerSide, filaments);
    }
    size_t N = filaments.size();
    if (N == 0) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "WindingLossesPEEC: no filaments generated");
    }

    // --- resistivity per filament (per unit length R = rho / area) ---
    auto resistivityModel = ResistivityModel::factory(ResistivityModels::WIRE_MATERIAL);
    std::vector<double> R(N);
    for (size_t k = 0; k < N; ++k) {
        auto wire = coil.resolve_wire(filaments[k].windingIndex);
        double rho = (*resistivityModel).get_resistivity(wire.resolve_material(), temperature);
        R[k] = rho / filaments[k].area; // [Ohm/m]
    }

    // --- ferrite-core boundary condition via CoilMesher's method of images ---
    // For each filament, get its image set (positions + (mu-n)/(mu+n) weights)
    // by reusing the exact rectangular-window image expansion already used by
    // the analytical field models. The REAL conductor (the (0,0) reflection,
    // which coincides with the filament centroid) is EXCLUDED here and handled
    // separately by the quadrature kernel below — the mesher returns it in the
    // middle of the list, not first, so it must be matched by position.
    // PEEC needs the core images for the 2D partial-inductance system to be
    // well-posed (without them the constrained log-kernel solve is singular and
    // collapses). Force at least one mirror level regardless of the caller's
    // field-mirroring setting; restore it afterwards.
    int savedMirroring = settings.get_magnetic_field_mirroring_dimension();
    if (savedMirroring < 1) {
        settings.set_magnetic_field_mirroring_dimension(1);
    }
    CoilMesherCenterModel imager;
    std::vector<std::vector<std::pair<double, double>>> imagePos(N); // (x,y)
    std::vector<std::vector<double>> imageWeight(N);
    {
        auto firstWire = coil.resolve_wire(filaments[0].windingIndex);
        for (size_t k = 0; k < N; ++k) {
            Turn synth;
            synth.set_coordinates({filaments[k].x, filaments[k].y});
            auto pts = imager.generate_mesh_inducing_turn(synth, firstWire, k, filaments[k].turnLength, core);
            for (auto& p : pts) {
                double px = p.get_point()[0], py = p.get_point()[1];
                if (std::abs(px - filaments[k].x) < 1e-12 && std::abs(py - filaments[k].y) < 1e-12) {
                    continue; // real conductor, handled by the quadrature term
                }
                imagePos[k].push_back({px, py});
                imageWeight[k].push_back(p.get_value());
            }
        }
    }
    if (savedMirroring < 1) {
        settings.set_magnetic_field_mirroring_dimension(savedMirroring);
    }

    // reference length to keep ln() dimensionless (cancels in the constrained solve)
    double refLength = core.get_winding_window().get_height().value();
    if (refLength <= 0) {
        refLength = 1.0;
    }

    // --- partial inductance matrix L [H/m] (real, frequency-independent) ---
    // L_ij = -(mu0/2pi) * [ <ln d>_ij(real, weight 1) + sum_images w*ln(d/ref) ]
    // Real term: 2D quadrature integral (i!=j) or self-GMD (i==j) -> convergent.
    // Image terms: far (window-scale), so the centroid point-log is accurate.
    double lnRef = std::log(refLength);
    Eigen::MatrixXd L(N, N);
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            double realLn;
            if (i == j) {
                realLn = std::log(self_gmd(filaments[i].width, filaments[i].height)) - lnRef;
            }
            else {
                realLn = avg_ln_between_cells(filaments[i].x, filaments[i].y, filaments[i].width, filaments[i].height,
                                              filaments[j].x, filaments[j].y, filaments[j].width, filaments[j].height) - lnRef;
            }
            double acc = realLn; // real conductor weight is 1
            for (size_t m = 0; m < imagePos[j].size(); ++m) {
                double dx = filaments[i].x - imagePos[j][m].first;
                double dy = filaments[i].y - imagePos[j][m].second;
                double d = std::max(std::sqrt(dx * dx + dy * dy), 1e-12);
                acc += imageWeight[j][m] * std::log(d / refLength);
            }
            L(i, j) = -(MU0 / (2 * std::numbers::pi)) * acc;
        }
    }

    // turn -> winding map and a current-direction sign per winding (primary +,
    // others -, matching the analytical convention).
    size_t numberWindings = coil.get_functional_description().size();
    std::vector<int> windingSign(numberWindings, -1);
    if (numberWindings > 0) {
        windingSign[0] = 1;
    }

    // incidence: filament -> turn (each turn is a parallel group of filaments)
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(N, numberTurns);
    for (size_t k = 0; k < N; ++k) {
        P(k, filaments[k].turnIndex) = 1.0;
    }

    // --- harmonic loop ---
    auto excitations = operatingPoint.get_excitations_per_winding();
    // gather a common harmonic list from winding 0's current
    auto refCurrent = excitations[0].get_current();
    if (!refCurrent || !refCurrent->get_harmonics()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "WindingLossesPEEC: operating point has no current harmonics");
    }
    auto refHarmonics = refCurrent->get_harmonics().value();
    auto frequencies = refHarmonics.get_frequencies();
    size_t numberHarmonics = frequencies.size();

    // amplitude reference for thresholding (fundamental of winding 0)
    double maxAmplitude = 0;
    for (auto& exc : excitations) {
        if (exc.get_current() && exc.get_current()->get_harmonics()) {
            for (double a : exc.get_current()->get_harmonics()->get_amplitudes()) {
                maxAmplitude = std::max(maxAmplitude, std::fabs(a));
            }
        }
    }

    double totalLosses = 0;
    std::vector<double> lossesPerWinding(numberWindings, 0.0);

    size_t systemSize = N + numberTurns;

    for (size_t h = 0; h < numberHarmonics; ++h) {
        double frequency = frequencies[h];

        // per-turn target current for this harmonic (signed by winding)
        Eigen::VectorXcd Iturns = Eigen::VectorXcd::Zero(numberTurns);
        double thisHarmonicMaxAmp = 0;
        for (size_t t = 0; t < numberTurns; ++t) {
            size_t w = coil.get_winding_index_by_name(turns[t].get_winding());
            auto exc = excitations[w];
            if (!exc.get_current() || !exc.get_current()->get_harmonics()) {
                continue;
            }
            auto amps = exc.get_current()->get_harmonics()->get_amplitudes();
            if (h >= amps.size()) {
                continue;
            }
            double amp = amps[h];
            // h==0 is DC; AC harmonics use rms = amplitude/sqrt(2)
            double current = (h == 0) ? amp : amp / std::sqrt(2.0);
            thisHarmonicMaxAmp = std::max(thisHarmonicMaxAmp, std::fabs(amp));
            Iturns(t) = cd(current * windingSign[w], 0.0);
        }

        if (frequency > 0 && maxAmplitude > 0 && thisHarmonicMaxAmp / maxAmplitude < _harmonicAmplitudeThreshold) {
            continue; // negligible harmonic
        }

        double omega = 2 * std::numbers::pi * frequency;

        // assemble [[Z, -P],[P^T, 0]] [I;V] = [0; Iturns], Z = R + jwL
        Eigen::MatrixXcd A = Eigen::MatrixXcd::Zero(systemSize, systemSize);
        for (size_t i = 0; i < N; ++i) {
            A(i, i) += cd(R[i], 0.0);
            for (size_t j = 0; j < N; ++j) {
                A(i, j) += cd(0.0, omega) * cd(L(i, j), 0.0);
            }
        }
        for (size_t i = 0; i < N; ++i) {
            for (size_t t = 0; t < numberTurns; ++t) {
                A(i, N + t) = cd(-P(i, t), 0.0);
                A(N + t, i) = cd(P(i, t), 0.0);
            }
        }
        Eigen::VectorXcd rhs = Eigen::VectorXcd::Zero(systemSize);
        for (size_t t = 0; t < numberTurns; ++t) {
            rhs(N + t) = Iturns(t);
        }

        // ColPivHouseholderQR is more robust than LU for the (mildly
        // ill-conditioned) 2D log-kernel PEEC system.
        Eigen::VectorXcd sol = A.colPivHouseholderQr().solve(rhs);

        // loss this harmonic = sum_k R_k |I_k|^2 * turnLength_k
        for (size_t k = 0; k < N; ++k) {
            cd Ik = sol(k);
            double p = R[k] * std::norm(Ik) * filaments[k].turnLength;
            totalLosses += p;
            lossesPerWinding[filaments[k].windingIndex] += p;
        }
    }

    // --- assemble output ---
    WindingLossesOutput out;
    out.set_method_used("PEEC_2D");
    out.set_origin(ResultOrigin::SIMULATION);
    out.set_winding_losses(totalLosses);
    out.set_temperature(temperature);

    std::vector<WindingLossesPerElement> perWinding;
    for (size_t w = 0; w < numberWindings; ++w) {
        WindingLossesPerElement el;
        OhmicLosses ohmic; // PEEC lumps skin+prox+dc together; report under ohmic for now
        ohmic.set_losses(lossesPerWinding[w]);
        ohmic.set_method_used("PEEC_2D");
        ohmic.set_origin(ResultOrigin::SIMULATION);
        el.set_ohmic_losses(ohmic);
        perWinding.push_back(el);
    }
    out.set_winding_losses_per_winding(perWinding);

    return out;
}

} // namespace OpenMagnetics
