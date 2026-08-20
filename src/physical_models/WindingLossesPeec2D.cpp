#include "physical_models/WindingLossesPeec2D.h"
#include "physical_models/WindingLosses.h"
#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/Reluctance.h"
#include "physical_models/Resistivity.h"
#include "support/CoilMesher.h"
#include "support/Settings.h"
#include "support/Exceptions.h"
#include "support/Utils.h"

#include <Eigen/Dense>
#include <complex>
#include <numbers>

namespace OpenMagnetics {

namespace {

constexpr double kMu0 = 4e-7 * std::numbers::pi;

// One rectangular sub-conductor of one turn's cross-section (per-metre 2D world).
struct Cell {
    double x = 0;          // centre
    double y = 0;
    double w = 0;          // dimensions
    double h = 0;
    double area = 0;       // conducting area this cell is billed for (round wires: masked share)
    size_t conductorIndex = 0;
};

// -----------------------------------------------------------------------------
// ln(geometric mean distance) machinery.
//
// Self: Rosa's classical approximation for a rectangle w x h,
//   GMD_self = 0.2235 * (w + h)
// (error < 0.3% over the aspect ratios the mesher produces; the end-to-end
// validation against the exact Bessel skin solution is the real acceptance).
//
// Mutual: exact only in the far field. ln r is smooth once the separation
// exceeds the cell sizes, so:
//   - far pairs (centre distance > 4x the summed half-diagonals): ln(centre
//     distance) — the multipole corrections are O((size/d)^2) there;
//   - near pairs: 4x4 Gauss-Legendre product quadrature on EACH rectangle
//     (256 point pairs) of the exact definition ln g = <ln|r1 - r2|>. This is
//     built ONCE per geometry (the matrix is frequency-independent), so the
//     cost sits far off the hot path.
// -----------------------------------------------------------------------------

inline double lnGmdSelf(double w, double h) {
    return std::log(0.2235 * (w + h));
}

// 4-point Gauss-Legendre abscissae/weights on [-1/2, 1/2].
constexpr std::array<double, 4> kGaussX = {-0.4305681557970263, -0.16999052179242816,
                                            0.16999052179242816, 0.4305681557970263};
constexpr std::array<double, 4> kGaussW = {0.17392742256872693, 0.32607257743127305,
                                           0.32607257743127305, 0.17392742256872693};

double lnGmdMutualQuadrature(const Cell& a, double bx, double by, double bw, double bh) {
    double acc = 0;
    for (size_t i = 0; i < 4; ++i) {
        double x1 = a.x + kGaussX[i] * a.w;
        double wx1 = kGaussW[i];
        for (size_t j = 0; j < 4; ++j) {
            double y1 = a.y + kGaussX[j] * a.h;
            double wy1 = kGaussW[j];
            for (size_t k = 0; k < 4; ++k) {
                double x2 = bx + kGaussX[k] * bw;
                double wx2 = kGaussW[k];
                for (size_t l = 0; l < 4; ++l) {
                    double y2 = by + kGaussX[l] * bh;
                    double wy2 = kGaussW[l];
                    double r = std::hypot(x1 - x2, y1 - y2);
                    // Two distinct cells can still share a corner point; a zero
                    // distance at a single quadrature node would poison the sum.
                    r = std::max(r, 1e-12);
                    acc += wx1 * wy1 * wx2 * wy2 * std::log(r);
                }
            }
        }
    }
    return acc;
}

inline double lnGmdMutual(const Cell& a, double bx, double by, double bw, double bh) {
    double d = std::hypot(a.x - bx, a.y - by);
    double reach = 2.0 * (std::hypot(a.w, a.h) + std::hypot(bw, bh));
    if (d > reach) {
        return std::log(d);
    }
    return lnGmdMutualQuadrature(a, bx, by, bw, bh);
}

// -----------------------------------------------------------------------------
// The winding-window image lattice — SAME frame selection, SAME reflection
// formulas and SAME (mu - k)/(mu + k) weights as CoilMesherCenterModel::
// generate_mesh_inducing_turn (CoilMesher.cpp), so this engine and the
// analytical field path can never disagree about the core geometry.
// TODO(#836): factor the lattice out of the mesher into one shared helper
// before AUTO dispatch is enabled; kept local for the validation shot.
// -----------------------------------------------------------------------------
struct ImageTerm {
    double x, y;
    double weight;  // includes the (mu-k)/(mu+k) ring weight; the m=n=0 original is weight 1
};

struct WindowFrame {
    double leftEdgeX, bottomY, A, B;
};

WindowFrame select_window_frame(const Core& core, double turnX) {
    auto processedDescription = core.get_processed_description().value();
    auto windingWindows = processedDescription.get_winding_windows();
    WindingWindowElement windingWindow = windingWindows[0];
    if (windingWindows.size() > 1) {
        double bestDistance = std::numeric_limits<double>::max();
        for (auto& candidate : windingWindows) {
            if (!candidate.get_coordinates() || !candidate.get_width()) {
                continue;
            }
            double distance = std::abs(turnX - candidate.get_coordinates().value()[0]);
            if (distance < bestDistance) {
                bestDistance = distance;
                windingWindow = candidate;
            }
        }
    }
    WindowFrame frame;
    frame.A = windingWindow.get_width().value();
    frame.B = windingWindow.get_height().value();
    if (windingWindows.size() > 1 && windingWindow.get_coordinates()) {
        frame.leftEdgeX = windingWindow.get_coordinates().value()[0] - frame.A / 2;
        frame.bottomY = windingWindow.get_coordinates().value()[1] - frame.B / 2;
    }
    else {
        frame.leftEdgeX = core.get_columns()[0].get_width() / 2;
        frame.bottomY = -frame.B / 2;
    }
    return frame;
}

std::vector<ImageTerm> image_lattice(double x, double y, const WindowFrame& frame,
                                     int mirroringDimension, double corePermeability) {
    std::vector<ImageTerm> terms;
    const int M = mirroringDimension;
    const int N = mirroringDimension;
    double crossingA = x - frame.leftEdgeX;
    double crossingB = y - frame.bottomY;
    for (int m = -M; m <= M; ++m) {
        for (int n = -N; n <= N; ++n) {
            double weight = (corePermeability - std::max(fabs(m), fabs(n))) /
                            (corePermeability + std::max(fabs(m), fabs(n)));
            double a = (m % 2 == 0) ? m * frame.A + crossingA : m * frame.A + frame.A - crossingA;
            double b = (n % 2 == 0) ? n * frame.B + crossingB : n * frame.B + frame.B - crossingB;
            terms.push_back({a + frame.leftEdgeX, b + frame.bottomY, weight});
        }
    }
    return terms;
}

// Mutual coefficient between cell i and conductor rectangle j INCLUDING j's
// images: -(mu0/2pi) * sum_k w_k * lnGMD(i, image_k(j)). The additive constant
// of the 2D limit is deliberately absent (see the header).
double couplingCoefficient(const Cell& receiver, double srcX, double srcY, double srcW, double srcH,
                           const WindowFrame& frame, int mirroringDimension, double corePermeability) {
    double sum = 0;
    for (const auto& image : image_lattice(srcX, srcY, frame, mirroringDimension, corePermeability)) {
        sum += image.weight * lnGmdMutual(receiver, image.x, image.y, srcW, srcH);
    }
    return -(kMu0 / (2 * std::numbers::pi)) * sum;
}

} // namespace

double WindingLossesPeec2D::estimate_angular_coverage(const Core& core) {
    auto family = core.get_shape_family();
    // A true pot core is a closed shell: the winding is covered over the whole revolution
    // (bar the two small lead slots), which is what OMFEM's ray march reports (1.000 on
    // P 3.3/2.6). MKF's column model cannot express "annular shell" — it records the same
    // bounding-box lateral a PQ has — so the family is the honest discriminator here.
    if (family == CoreShapeFamily::P || family == CoreShapeFamily::PM) {
        return 1.0;
    }
    double covered = 0;
    // get_columns() returns BY VALUE — bind a named local before iterating, or the
    // range-for holds a dangling reference to a destroyed temporary (repo memory:
    // mas-optional-value-dangling-ref; it corrupted the heap here, not merely read
    // garbage).
    const auto columns = core.get_columns();
    for (const auto& column : columns) {
        if (column.get_type() != ColumnType::LATERAL) {
            continue;
        }
        double xCentre = std::abs(column.get_coordinates()[0]);
        double innerFace = std::max(1e-9, xCentre - column.get_width() / 2);
        double tangentialHalf = column.get_depth() / 2;
        covered += 2 * std::atan2(tangentialHalf, innerFace);
    }
    if (covered <= 0) {
        return 1.0;  // no lateral columns recorded: treat as fully covered rather than guess
    }
    return std::clamp(covered / (2 * std::numbers::pi), 0.0, 1.0);
}

WindingLossesOutput WindingLossesPeec2D::calculate_losses(Magnetic magnetic, OperatingPoint operatingPoint,
                                                          double temperature) {
    auto& settings = Settings::GetInstance();
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();

    if (!core.get_processed_description()) {
        core.process_data();
        magnetic.set_core(core);
    }
    if (core.get_shape_family() == CoreShapeFamily::T) {
        throw NotImplementedException(
            "Peec2D winding losses are only implemented for concentric windows: a toroidal (round) "
            "window has no rectangular image frame. Use the analytical winding-loss path for toroids.");
    }
    if (!coil.get_turns_description()) {
        throw CoilNotProcessedException("Peec2D winding losses need a wound coil (turns description)");
    }

    // Normalise the operating point exactly like WindingLosses::calculate_losses does.
    bool needsProcessing = false;
    for (const auto& excitation : operatingPoint.get_excitations_per_winding()) {
        if (excitation.get_current() &&
            (!excitation.get_current()->get_waveform() || !excitation.get_current()->get_harmonics() ||
             !excitation.get_current()->get_processed() ||
             !excitation.get_current()->get_processed()->get_rms())) {
            needsProcessing = true;
            break;
        }
    }
    if (needsProcessing) {
        auto magnetizingInductance = resolve_dimensional_values(
            MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(
                magnetic.get_core(), magnetic.get_coil()).get_magnetizing_inductance());
        operatingPoint = Inputs::process_operating_point(operatingPoint, magnetizingInductance,
                                                         magnetic.get_mutable_coil().get_turns_ratios());
    }

    // Stage 1, unchanged: DC/ohmic losses (also owns the connection copper).
    auto windingLossesOutput = WindingOhmicLosses::calculate_ohmic_losses(coil, operatingPoint, temperature);

    const auto turns = coil.get_turns_description().value();
    const auto wirePerWinding = coil.get_wires();

    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        if (coil.get_functional_description()[windingIndex].get_number_parallels() != 1) {
            throw NotImplementedException(
                "Peec2D winding losses: winding '" + coil.get_functional_description()[windingIndex].get_name() +
                "' has " + std::to_string(coil.get_functional_description()[windingIndex].get_number_parallels()) +
                " parallels. Solving the parallel current split is the designed extension of the "
                "constraint machinery but is not implemented yet; use the analytical path.");
        }
        if (coil.get_wire_type(windingIndex) == WireType::LITZ) {
            throw NotImplementedException(
                "Peec2D winding losses: litz keeps the analytical path (a solid-equivalent PEEC "
                "would misrepresent strand-level behaviour).");
        }
    }

    // ---- Harmonic selection: the SAME pruning the analytical pipeline uses. ----
    auto commonHarmonicIndexes =
        CoilMesher().get_common_harmonic_indexes(operatingPoint, settings.get_harmonic_amplitude_threshold());
    auto primaryHarmonics =
        operatingPoint.get_excitations_per_winding()[0].get_current()->get_harmonics().value();

    // Winding current direction convention: primary positive, the rest negative
    // (same default as CoilMesher::generate_mesh_inducing_coil).
    std::vector<double> currentDirectionPerWinding(coil.get_functional_description().size(), -1.0);
    currentDirectionPerWinding[0] = 1.0;

    // Per-winding harmonic amplitudes on the COMMON harmonic grid.
    std::vector<Harmonics> harmonicsPerWinding;
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        harmonicsPerWinding.push_back(operatingPoint.get_excitations_per_winding()[windingIndex]
                                          .get_current()->get_harmonics().value());
    }

    // ---- Resistivity and skin depth bounds. ----
    auto resistivityModel = ResistivityModel::factory(ResistivityModels::WIRE_MATERIAL);
    double highestFrequency = 0;
    for (auto harmonicIndex : commonHarmonicIndexes) {
        double frequency = primaryHarmonics.get_frequencies()[harmonicIndex];
        if (frequency > 0) {
            highestFrequency = std::max(highestFrequency, frequency);
        }
    }
    if (highestFrequency <= 0) {
        // Pure DC excitation: the ohmic stage already covers it.
        return windingLossesOutput;
    }

    // ---- Gap conductors: per FUNCTIONAL gap, MMF(h) = phi(h) * R_gap. ----
    // phi(h) = L * I_mag(h) / N with the same inductance model the pipeline uses.
    struct GapConductor {
        double x, y;
        double length;                    // the physical gap length: the MMF source's extent
        double mmfPerUnitPrimaryCurrent;  // A of MMF per A of primary harmonic current
    };
    std::vector<GapConductor> gapConductors;
    {
        auto gapping = core.get_functional_description().get_gapping();
        bool anyFunctional = false;
        for (auto& gap : gapping) {
            if (gap.get_type() == GapType::SUBTRACTIVE || gap.get_type() == GapType::ADDITIVE) {
                anyFunctional = true;
            }
        }
        if (anyFunctional) {
            double numberPrimaryTurns =
                double(coil.get_functional_description()[0].get_number_turns());
            double magnetizingInductance = resolve_dimensional_values(
                MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(
                    core, coil).get_magnetizing_inductance());
            auto reluctanceModel = ReluctanceModel::factory();
            for (auto& gap : gapping) {
                if (gap.get_type() != GapType::SUBTRACTIVE && gap.get_type() != GapType::ADDITIVE) {
                    continue;  // residual mating surfaces do not fringe (ABT #832)
                }
                if (!gap.get_coordinates() || gap.get_coordinates().value()[0] < 0) {
                    continue;  // same x>=0 convention as the fringing path
                }
                double gapReluctance = reluctanceModel->get_gap_reluctance(gap).get_reluctance();
                // flux per ampere of primary current: L * 1 / N
                double fluxPerAmp = magnetizingInductance / numberPrimaryTurns;
                GapConductor gapConductor;
                // The fictitious conductor lives ON the winding-window boundary, at the
                // gap's surface ADJACENT to the window — not at the gap's own centre
                // coordinates. A centre-column gap is recorded at x = 0 (the column
                // AXIS), which lies outside the 2D window frame entirely: placed there
                // it couples ~nothing and its images wrap around a frame it is not in
                // (first validation run: foil crowding read 0.10-0.15x of FEM). The
                // Roshen fringing model makes the same adjacency choice (distance from
                // the gap CORNER, gap x +- sectionDimensions/2).
                double gapX = gap.get_coordinates().value()[0];
                double halfSection = gap.get_section_dimensions()
                    ? gap.get_section_dimensions().value()[0] / 2 : 0.0;
                GapConductor placed;
                placed.length = gap.get_length();
                placed.y = gap.get_coordinates().value()[1];
                if (std::abs(gapX) < halfSection) {
                    // centre column: the window-adjacent surface is the column edge
                    placed.x = gapX + halfSection;
                }
                else {
                    // lateral column at +x: the window sits to its LEFT
                    placed.x = gapX - halfSection;
                }
                placed.mmfPerUnitPrimaryCurrent = fluxPerAmp * gapReluctance;
                gapConductor = placed;
                gapConductors.push_back(gapConductor);
            }
        }
    }

    // ---- Cell mesh. ----
    std::vector<Cell> cells;
    std::vector<double> conductorLength(turns.size(), 0.0);
    std::vector<size_t> conductorWinding(turns.size(), 0);
    std::vector<size_t> cellsBeginPerConductor(turns.size() + 1, 0);

    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        const auto& turn = turns[turnIndex];
        size_t windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        conductorWinding[turnIndex] = windingIndex;
        conductorLength[turnIndex] = turn.get_length();
        auto wire = coil.resolve_wire(windingIndex);
        double resistivity = (*resistivityModel).get_resistivity(wire.resolve_material(), temperature);
        (void)resistivity;

        double turnX = turn.get_coordinates()[0];
        double turnY = turn.get_coordinates()[1];
        cellsBeginPerConductor[turnIndex] = cells.size();

        auto pushGrid = [&](double width, double height, bool maskRound, double diameter) {
            double deltaMin = std::sqrt(
                (*resistivityModel).get_resistivity(wire.resolve_material(), temperature) /
                (std::numbers::pi * highestFrequency * kMu0));
            double thin = std::min(width, height);
            double wide = std::max(width, height);
            // BOTH dimensions resolve the skin depth: current crowds along the WIDE
            // dimension's edges too (the 2D bar problem), and under-resolving it is
            // indistinguishable from missing physics — the first validation run read
            // 25% low on an isolated 4:1 bar with 2-skin-depth-wide cells.
            size_t nThin = std::clamp<size_t>(size_t(std::ceil(cellsPerSkinDepth * thin / deltaMin)),
                                              minimumCellsThin, maximumCellsThin);
            size_t nWide = std::clamp<size_t>(size_t(std::ceil(cellsPerSkinDepth * wide / deltaMin)),
                                              minimumCellsWide, maximumCellsWide);
            // For round wires the grid is square-ish: both dimensions penetrate.
            size_t nx, ny;
            if (maskRound) {
                size_t n = std::clamp<size_t>(size_t(std::ceil(cellsPerSkinDepth * diameter / deltaMin)),
                                              6, 24);
                nx = n;
                ny = n;
                width = diameter;
                height = diameter;
            }
            else if (width >= height) {
                nx = nWide;
                ny = nThin;
            }
            else {
                nx = nThin;
                ny = nWide;
            }
            // Refinement near a functional gap (the paper's 40x9-near-gap pattern):
            // the crowding profile lives along the WIDE dimension at the gap's height,
            // so it is the wide-axis cell count that doubles — the first version
            // doubled nx unconditionally, which for a vertical foil is the 0.1 mm
            // THICKNESS axis, and used a trigger distance the old (mis-placed) gap
            // coordinates could never satisfy.
            for (const auto& gapConductor : gapConductors) {
                double trigger = 3e-3 + std::max(width, height) / 2;
                if (std::hypot(turnX - gapConductor.x, turnY - gapConductor.y) < trigger) {
                    if (maskRound) {
                        nx = std::min<size_t>(nx * 2, 32);
                        ny = std::min<size_t>(ny * 2, 32);
                    }
                    else if (width >= height) {
                        nx = std::min<size_t>(nx * 2, maximumCellsWide * 2);
                    }
                    else {
                        ny = std::min<size_t>(ny * 2, maximumCellsWide * 2);
                    }
                    break;
                }
            }
            double dx = width / double(nx);
            double dy = height / double(ny);
            double keptArea = 0;
            size_t firstCell = cells.size();
            for (size_t i = 0; i < nx; ++i) {
                for (size_t j = 0; j < ny; ++j) {
                    double cx = turnX - width / 2 + (i + 0.5) * dx;
                    double cy = turnY - height / 2 + (j + 0.5) * dy;
                    if (maskRound &&
                        std::hypot(cx - turnX, cy - turnY) > diameter / 2) {
                        continue;
                    }
                    Cell cell;
                    cell.x = cx;
                    cell.y = cy;
                    cell.w = dx;
                    cell.h = dy;
                    cell.area = dx * dy;
                    cell.conductorIndex = turnIndex;
                    cells.push_back(cell);
                    keptArea += cell.area;
                }
            }
            // Preserve the exact conducting area (masked circles): scale the cell
            // areas so the DC resistance of the mesh equals the wire's.
            double trueArea = maskRound ? std::numbers::pi * diameter * diameter / 4 : width * height;
            double areaScale = trueArea / keptArea;
            for (size_t c = firstCell; c < cells.size(); ++c) {
                cells[c].area *= areaScale;
            }
        };

        switch (wire.get_type()) {
            case WireType::ROUND: {
                double d = resolve_dimensional_values(wire.get_conducting_diameter().value());
                pushGrid(d, d, true, d);
                break;
            }
            case WireType::RECTANGULAR:
            case WireType::FOIL:
            case WireType::PLANAR: {
                double w = resolve_dimensional_values(wire.get_conducting_width().value());
                double h = resolve_dimensional_values(wire.get_conducting_height().value());
                pushGrid(w, h, false, 0);
                break;
            }
            default:
                throw NotImplementedException("Peec2D winding losses: unsupported wire type");
        }
    }
    cellsBeginPerConductor[turns.size()] = cells.size();

    const size_t numberCells = cells.size();
    const size_t numberConductors = turns.size();
    if (numberCells > maximumCells) {
        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION,
            "Peec2D winding losses: the mesh needs " + std::to_string(numberCells) +
            " cells, above the configured cap of " + std::to_string(maximumCells) +
            ". Raise WindingLossesPeec2D::maximumCells or use the analytical path — a silently "
            "coarsened mesh would misreport losses.");
    }

    _diagnostics = Diagnostics{};
    _diagnostics.totalCells = numberCells;
    _diagnostics.totalConductors = numberConductors;
    _diagnostics.perTurnPerHarmonic.resize(numberConductors);

    // ---- Frequency-independent matrices. ----
    const int mirroringDimension = settings.get_magnetic_field_mirroring_dimension();
    const double corePermeability = core.get_initial_permeability(temperature);
    WindowFrame frame = select_window_frame(core, turns[0].get_coordinates()[0]);

    // ABT #837: two inductance matrices — the fully-imaged window and the open-air
    // (image-free) one. See the angularCoverage comment in the header: the losses of the
    // two solves are blended by the core's angular coverage, because an E/PQ/ETD winding
    // is only backed by core over a fraction of its revolution.
    const double coverage = angularCoverage ? std::clamp(angularCoverage.value(), 0.0, 1.0) : 1.0;
    const bool needsOpenSolve = coverage < 0.999;

    Eigen::MatrixXd inductanceMatrix(numberCells, numberCells);
    Eigen::MatrixXd inductanceMatrixOpen;
    if (needsOpenSolve) {
        inductanceMatrixOpen.resize(numberCells, numberCells);
    }
    for (size_t i = 0; i < numberCells; ++i) {
        for (size_t j = i; j < numberCells; ++j) {
            double value = couplingCoefficient(cells[i], cells[j].x, cells[j].y, cells[j].w, cells[j].h,
                                               frame, mirroringDimension, corePermeability);
            inductanceMatrix(i, j) = value;
            inductanceMatrix(j, i) = value;
            if (needsOpenSolve) {
                double bare = (i == j)
                    ? -(kMu0 / (2 * std::numbers::pi)) * lnGmdSelf(cells[i].w, cells[i].h)
                    : -(kMu0 / (2 * std::numbers::pi)) * lnGmdMutual(cells[i], cells[j].x, cells[j].y,
                                                                     cells[j].w, cells[j].h);
                inductanceMatrixOpen(i, j) = bare;
                inductanceMatrixOpen(j, i) = bare;
            }
        }
    }
    // Gap coupling column per gap: -(mu0/2pi) sum_images w * ln r(cell, gap images).
    Eigen::MatrixXd gapCoupling(numberCells, gapConductors.size());
    Eigen::MatrixXd gapCouplingOpen(numberCells, gapConductors.size());
    for (size_t g = 0; g < gapConductors.size(); ++g) {
        for (size_t i = 0; i < numberCells; ++i) {
            if (needsOpenSolve) {
                Cell gapCell;
                gapCell.x = gapConductors[g].x; gapCell.y = gapConductors[g].y;
                gapCell.w = gapConductors[g].length / 2; gapCell.h = gapConductors[g].length;
                gapCouplingOpen(i, g) = -(kMu0 / (2 * std::numbers::pi)) *
                    lnGmdMutual(cells[i], gapCell.x, gapCell.y, gapCell.w, gapCell.h);
            }
            // The MMF source carries the PHYSICAL extent of the gap, not a filament.
            // Kovacevic's replacement assumes the gap is small against the winding
            // distance; a foil 1.6 mm from a 1 mm gap violates that, and a 1 um
            // filament there over-drove the crowding 1.5x at 50 kHz growing to 3.9x
            // at 500 kHz. A gap-length-sized source rectangle bounds the near field
            // the way the real gap aperture does (the GMD quadrature spreads it).
            gapCoupling(i, g) = couplingCoefficient(cells[i], gapConductors[g].x, gapConductors[g].y,
                                                    gapConductors[g].length / 2, gapConductors[g].length,
                                                    frame, mirroringDimension, corePermeability);
        }
    }

    std::vector<double> cellResistancePerMetre(numberCells);
    for (size_t i = 0; i < numberCells; ++i) {
        auto wire = coil.resolve_wire(conductorWinding[cells[i].conductorIndex]);
        double resistivity = (*resistivityModel).get_resistivity(wire.resolve_material(), temperature);
        cellResistancePerMetre[i] = resistivity / cells[i].area;
    }

    // Per-conductor DC resistance per metre (for the dc-equivalent split).
    std::vector<double> conductorDcResistancePerMetre(numberConductors);
    for (size_t c = 0; c < numberConductors; ++c) {
        double conductance = 0;
        for (size_t i = cellsBeginPerConductor[c]; i < cellsBeginPerConductor[c + 1]; ++i) {
            conductance += 1.0 / cellResistancePerMetre[i];
        }
        conductorDcResistancePerMetre[c] = 1.0 / conductance;
    }

    // ---- Isolated-turn factorisations (block solves), cached per distinct geometry. ----
    // Isolated = the turn's own cells only, free space (no images, no gap): the
    // classical skin-effect problem, solved with the same machinery.
    auto solveIsolated = [&](size_t conductorIdx, double omega) -> double {
        size_t begin = cellsBeginPerConductor[conductorIdx];
        size_t end = cellsBeginPerConductor[conductorIdx + 1];
        size_t n = end - begin;
        Eigen::MatrixXcd system(n + 1, n + 1);
        system.setZero();
        for (size_t i = 0; i < n; ++i) {
            const Cell& cellI = cells[begin + i];
            for (size_t j = 0; j < n; ++j) {
                const Cell& cellJ = cells[begin + j];
                double inductance = (i == j)
                    ? -(kMu0 / (2 * std::numbers::pi)) * lnGmdSelf(cellI.w, cellI.h)
                    : -(kMu0 / (2 * std::numbers::pi)) * lnGmdMutual(cellI, cellJ.x, cellJ.y, cellJ.w, cellJ.h);
                system(i, j) = std::complex<double>(0, omega * inductance);
            }
            system(i, i) += cellResistancePerMetre[begin + i];
            system(i, n) = 1.0;
            system(n, i) = 1.0;
        }
        Eigen::VectorXcd rhs = Eigen::VectorXcd::Zero(n + 1);
        rhs(n) = 1.0;  // unit peak current
        Eigen::VectorXcd solution = system.partialPivLu().solve(rhs);
        double lossPerMetrePerAmpSquared = 0;
        for (size_t i = 0; i < n; ++i) {
            lossPerMetrePerAmpSquared += 0.5 * cellResistancePerMetre[begin + i] * std::norm(solution(i));
        }
        return lossPerMetrePerAmpSquared;  // W/m at 1 A peak
    };

    // ---- Per-harmonic full solves. ----
    const size_t systemSize = numberCells + numberConductors;
    Eigen::MatrixXcd system(systemSize, systemSize);
    Eigen::VectorXcd rhs(systemSize);

    // Output accumulators, schema-shaped.
    auto windingLossesPerTurn = windingLossesOutput.get_winding_losses_per_turn().value();
    std::vector<LossElementPerHarmonic> skinPerTurn(numberConductors);
    std::vector<LossElementPerHarmonic> proximityPerTurn(numberConductors);
    for (size_t c = 0; c < numberConductors; ++c) {
        skinPerTurn[c].set_method_used("Peec2D");
        skinPerTurn[c].set_origin(ResultOrigin::SIMULATION);
        skinPerTurn[c].get_mutable_harmonic_frequencies().push_back(0);
        skinPerTurn[c].get_mutable_losses_per_harmonic().push_back(0);
        proximityPerTurn[c].set_method_used("Peec2D");
        proximityPerTurn[c].set_origin(ResultOrigin::SIMULATION);
        proximityPerTurn[c].get_mutable_harmonic_frequencies().push_back(0);
        proximityPerTurn[c].get_mutable_losses_per_harmonic().push_back(0);
    }
    double totalExtraLosses = 0;

    for (auto harmonicIndex : commonHarmonicIndexes) {
        double frequency = primaryHarmonics.get_frequencies()[harmonicIndex];
        if (frequency <= 0) {
            continue;  // DC is the ohmic stage's job
        }
        double omega = 2 * std::numbers::pi * frequency;

        double primaryAmplitude = harmonicsPerWinding[0].get_amplitudes()[harmonicIndex];

        // One bordered solve for a given coupling model; returns the per-conductor
        // per-metre losses.
        auto solveWith = [&](const Eigen::MatrixXd& L, const Eigen::MatrixXd& gapL) {
            system.setZero();
            for (size_t i = 0; i < numberCells; ++i) {
                for (size_t j = 0; j < numberCells; ++j) {
                    system(i, j) = std::complex<double>(0, omega * L(i, j));
                }
                system(i, i) += cellResistancePerMetre[i];
                system(i, numberCells + cells[i].conductorIndex) = 1.0;
            }
            for (size_t c = 0; c < numberConductors; ++c) {
                for (size_t i = cellsBeginPerConductor[c]; i < cellsBeginPerConductor[c + 1]; ++i) {
                    system(numberCells + c, i) = 1.0;
                }
            }
            rhs.setZero();
            for (size_t g = 0; g < gapConductors.size(); ++g) {
                double gapCurrent = gapConductors[g].mmfPerUnitPrimaryCurrent * primaryAmplitude;
                for (size_t i = 0; i < numberCells; ++i) {
                    rhs(i) -= std::complex<double>(0, omega * gapL(i, g) * gapCurrent);
                }
            }
            for (size_t c = 0; c < numberConductors; ++c) {
                size_t windingIndex = conductorWinding[c];
                double amplitude = harmonicsPerWinding[windingIndex].get_amplitudes()[harmonicIndex];
                rhs(numberCells + c) = currentDirectionPerWinding[windingIndex] * amplitude;
            }
            Eigen::VectorXcd solution = system.partialPivLu().solve(rhs);
            std::vector<double> lossPerConductor(numberConductors, 0.0);
            for (size_t c = 0; c < numberConductors; ++c) {
                for (size_t i = cellsBeginPerConductor[c]; i < cellsBeginPerConductor[c + 1]; ++i) {
                    lossPerConductor[c] += 0.5 * cellResistancePerMetre[i] * std::norm(solution(i));
                }
            }
            return lossPerConductor;
        };

        std::vector<double> coveredLoss = solveWith(inductanceMatrix, gapCoupling);
        std::vector<double> openLoss = needsOpenSolve ? solveWith(inductanceMatrixOpen, gapCouplingOpen)
                                                      : coveredLoss;

        // Losses.
        for (size_t c = 0; c < numberConductors; ++c) {
            // 2.5D angular blend (see the header): each azimuthal sector sees its own 2D
            // problem — backed by core over the covered fraction, open air elsewhere.
            double fullLossPerMetre = coverage * coveredLoss[c] + (1 - coverage) * openLoss[c];
            size_t windingIndex = conductorWinding[c];
            double amplitude = harmonicsPerWinding[windingIndex].get_amplitudes()[harmonicIndex];
            double dcEquivalentPerMetre = 0.5 * conductorDcResistancePerMetre[c] * amplitude * amplitude;
            double isolatedPerMetre = solveIsolated(c, omega) * amplitude * amplitude;

            double length = conductorLength[c];
            double fullLoss = fullLossPerMetre * length;
            double isolatedLoss = isolatedPerMetre * length;
            double dcEquivalentLoss = dcEquivalentPerMetre * length;

            // The classical decomposition, by definition. Numerical guards: the
            // full solve can undercut the isolated one only by round-off.
            double skinLoss = std::max(0.0, isolatedLoss - dcEquivalentLoss);
            double proximityLoss = fullLoss - isolatedLoss;

            skinPerTurn[c].get_mutable_harmonic_frequencies().push_back(frequency);
            skinPerTurn[c].get_mutable_losses_per_harmonic().push_back(skinLoss);
            proximityPerTurn[c].get_mutable_harmonic_frequencies().push_back(frequency);
            proximityPerTurn[c].get_mutable_losses_per_harmonic().push_back(proximityLoss);
            totalExtraLosses += skinLoss + proximityLoss;

            _diagnostics.perTurnPerHarmonic[c].push_back(
                TurnHarmonicLoss{frequency, fullLoss, isolatedLoss, dcEquivalentLoss});
        }
    }

    for (size_t c = 0; c < numberConductors; ++c) {
        windingLossesPerTurn[c].set_skin_effect_losses(skinPerTurn[c]);
        windingLossesPerTurn[c].set_proximity_effect_losses(proximityPerTurn[c]);
    }
    windingLossesOutput.set_winding_losses_per_turn(windingLossesPerTurn);
    windingLossesOutput.set_winding_losses(windingLossesOutput.get_winding_losses() + totalExtraLosses);
    windingLossesOutput.set_method_used("Peec2D");
    windingLossesOutput = WindingLosses::combine_turn_losses(windingLossesOutput, coil);
    return windingLossesOutput;
}

} // namespace OpenMagnetics
