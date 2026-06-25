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
#include <numbers>

namespace OpenMagnetics {

using cd = std::complex<double>;
static const double MU0 = 4 * std::numbers::pi * 1e-7;

// Geometric mean distance of a rectangle to itself (Rosa), used for the
// self partial inductance term so ln(distance) doesn't blow up at d=0.
static double self_gmd(double w, double h) {
    return 0.2235 * (w + h);
}

// Subdivide a single turn's conducting cross-section into filaments.
static void mesh_turn_into_filaments(const Turn& turn, Wire wire, size_t turnIndex, size_t windingIndex,
                                     size_t nPerDim, std::vector<WindingLossesPEEC::Filament>& out) {
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

    double dw = fullW / nPerDim;
    double dh = fullH / nPerDim;

    // First pass: collect candidate centroids (round wire keeps only cells whose
    // centre is inside the conductor) so the total filament area matches the
    // real conducting area.
    std::vector<std::pair<double, double>> centres;
    for (size_t i = 0; i < nPerDim; ++i) {
        for (size_t j = 0; j < nPerDim; ++j) {
            double fx = cx - fullW / 2 + (i + 0.5) * dw;
            double fy = cy - fullH / 2 + (j + 0.5) * dh;
            if (isRound) {
                double rx = fx - cx, ry = fy - cy;
                if (rx * rx + ry * ry > radius * radius) {
                    continue;
                }
            }
            centres.push_back({fx, fy});
        }
    }
    if (centres.empty()) {
        return;
    }

    double conductingArea = isRound ? (std::numbers::pi * radius * radius) : (fullW * fullH);
    double filamentArea = conductingArea / centres.size();

    for (auto& [fx, fy] : centres) {
        WindingLossesPEEC::Filament f;
        f.x = fx;
        f.y = fy;
        f.area = filamentArea;
        f.width = dw;
        f.height = dh;
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

    // --- mesh every turn into filaments; reject litz (impractical to mesh strands) ---
    std::vector<Filament> filaments;
    // Auto-reduce resolution so the dense solve stays tractable on big windings.
    size_t nPerDim = _filamentsPerDimension;
    while (nPerDim > 1 && numberTurns * nPerDim * nPerDim > 1300) {
        --nPerDim;
    }

    for (size_t turnIndex = 0; turnIndex < numberTurns; ++turnIndex) {
        auto turn = turns[turnIndex];
        size_t windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        auto wire = coil.resolve_wire(windingIndex);
        if (wire.get_type() == WireType::LITZ) {
            throw NotImplementedException("WindingLossesPEEC: litz must use the analytical model (strand meshing is impractical); homogenise the bundle to enable PEEC");
        }
        mesh_turn_into_filaments(turn, wire, turnIndex, windingIndex, nPerDim, filaments);
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
    // the analytical field models. The first entry (weight 1) is the filament
    // itself.
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
                imagePos[k].push_back({p.get_point()[0], p.get_point()[1]});
                imageWeight[k].push_back(p.get_value());
            }
        }
    }

    // reference length to keep ln() dimensionless (cancels in the constrained solve)
    double refLength = core.get_winding_window().get_height().value();
    if (refLength <= 0) {
        refLength = 1.0;
    }

    // --- partial inductance matrix L [H/m] (real, frequency-independent) ---
    // L_ij = -(mu0/2pi) * sum over images of j: w * ln(dist(i, image_of_j)/ref)
    Eigen::MatrixXd L(N, N);
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            double acc = 0;
            for (size_t m = 0; m < imagePos[j].size(); ++m) {
                double dx = filaments[i].x - imagePos[j][m].first;
                double dy = filaments[i].y - imagePos[j][m].second;
                double d = std::sqrt(dx * dx + dy * dy);
                bool isSelfReal = (i == j) && (m == 0); // real-on-real -> use GMD
                if (isSelfReal) {
                    d = self_gmd(filaments[j].width, filaments[j].height);
                }
                d = std::max(d, 1e-12);
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

        Eigen::VectorXcd sol = A.partialPivLu().solve(rhs);

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
