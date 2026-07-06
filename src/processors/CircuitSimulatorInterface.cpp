#include "support/Utils.h"
#include "Defaults.h"
#include "processors/CircuitSimulatorInterface.h"
#include "processors/CircuitSimulatorExporterHelpers.h"
#include "processors/Inputs.h"
#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/LeakageInductance.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingLosses.h"
#include "processors/Sweeper.h"
#include "support/Settings.h"
#include <filesystem>
#include <ctime>
// levmar.h removed - using Eigen LevenbergMarquardt
#include <float.h>
#include "support/Exceptions.h"
#include "physical_models/CoreLosses.h"
#include "physical_models/InitialPermeability.h"
#include "physical_models/StrayCapacitance.h"

namespace OpenMagnetics {

// Extract saturation parameters from a magnetic component
SaturationParameters get_saturation_parameters(const Magnetic& magnetic, double temperature) {
    SaturationParameters params;
    params.valid = false;
    params.Isat = 0;

    try {
        auto core = magnetic.get_core();
        auto coil = magnetic.get_coil();
        
        // Get saturation flux density
        params.Bsat = core.get_magnetic_flux_density_saturation(temperature, true);
        
        // Get initial permeability
        params.mu_r = core.get_initial_permeability(temperature);
        
        // Get effective parameters
        params.Ae = core.get_effective_area();
        params.le = core.get_effective_length();
        
        // Get primary turns
        params.primaryTurns = coil.get_functional_description()[0].get_number_turns();
        
        // Get magnetizing inductance
        params.Lmag = resolve_dimensional_values(
            MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(magnetic)
                .get_magnetizing_inductance());
        
        // Saturation current comes from the authoritative physical model
        // (Magnetic::calculate_saturation_current = gapped B_sat*N*A_e/L), NOT a
        // formula reimplemented here -- that is what drifted and produced the
        // ungapped, ~10x too-low I_sat in ABT #40. proportion=true matches the
        // proportioned B_sat used above and by the model internally. A mutable
        // copy is needed because calculate_saturation_current uses the mutable
        // core/coil accessors.
        if (params.Bsat > 0 && params.mu_r > 0 && params.Ae > 0 &&
            params.le > 0 && params.primaryTurns > 0 && params.Lmag > 0) {
            Magnetic mutableMagnetic = magnetic;
            params.Isat = mutableMagnetic.calculate_saturation_current(temperature, true);
            params.valid = params.Isat > 0;
        }
    } catch (const std::exception& e) {
        // Silently exporting a LINEAR inductor when saturation was requested hid
        // real data/model failures. Surface them instead.
        throw std::runtime_error(std::string("get_saturation_parameters failed (saturation modeling was requested): ") + e.what());
    }

    return params;
}

template<typename T> std::vector<T> convolution_valid(std::vector<T> const &f, std::vector<T> const &g) {
    int const nf = f.size();
    int const ng = g.size();
    std::vector<T> const &min_v = (nf < ng)? f : g;
    std::vector<T> const &max_v = (nf < ng)? g : f;
    int const n  = std::max(nf, ng) - std::min(nf, ng) + 1;
    std::vector<T> out(n, T());
    for(auto i(0); i < n; ++i) {
        for(int j(min_v.size() - 1), k(i); j >= 0; --j) {
            out[i] += min_v[j] * max_v[k];
            ++k;
        }
    }
    return out;
}

std::string to_string(double d, size_t precision) {
    std::ostringstream out;
    out.precision(precision);
    out << std::fixed << d;
    return out.str();
}

std::vector<std::string> to_string(std::vector<double> v, size_t precision=12) {
    std::vector<std::string> s;
    for (auto d : v){
        s.push_back(to_string(d, precision));
    }
    return s;
}

CircuitSimulatorExporter::CircuitSimulatorExporter(CircuitSimulatorExporterModels program) {
    _model = CircuitSimulatorExporterModel::factory(program);
}

CircuitSimulatorExporter::CircuitSimulatorExporter() {
    _model = CircuitSimulatorExporterModel::factory(CircuitSimulatorExporterModels::SIMBA);
}

double CircuitSimulatorExporter::analytical_model(double x[], double frequency) {
    return x[0] + x[1] * sqrt(frequency) + x[2] * frequency; 
}

void CircuitSimulatorExporter::analytical_func(double *p, double *x, int m, int n, void *data) {
    double* frequencies = static_cast <double*> (data);
    for(int i=0; i<n; ++i) {
        x[i]=analytical_model(p, frequencies[i]);
    }
}

std::complex<double> parallel(std::complex<double> z0, std::complex<double> z1) {
    return 1.0 / (1.0 / z0 + 1.0 / z1);
}

double CircuitSimulatorExporter::ladder_model(double x[], double frequency, double dcResistance) {
    double w = 2 * std::numbers::pi * frequency;
    auto R0 = std::complex<double>(dcResistance);
    auto R1 = std::complex<double>(x[0], 0);
    auto L1 = std::complex<double>(0, w * x[1]);
    auto R2 = std::complex<double>(x[2], 0);
    auto L2 = std::complex<double>(0, w * x[3]);
    auto R3 = std::complex<double>(x[4], 0);
    auto L3 = std::complex<double>(0, w * x[5]);
    auto R4 = std::complex<double>(x[6], 0);
    auto L4 = std::complex<double>(0, w * x[7]);
    auto R5 = std::complex<double>(x[8], 0);
    auto L5 = std::complex<double>(0, w * x[9]);

    // FIXED: Removed overly tight x[1] > 1e-7 constraint
    // FIXED: Return large penalty instead of 0 to guide optimizer
    for(int i=0; i<10; ++i) {
        if (x[i] < 0) {
            return 1e30;  // Large penalty value instead of 0
        }
    }
    // Reject ladder stages whose L/R SETTLING TIME is too long for a transient
    // sim to reach steady state. levmar can otherwise fit a physically-absurd
    // inductor (e.g. ~10 mH with a ~19 Ω stage -> ~550 µs settling) that shapes
    // the broadband curve on average but never settles in the transient regulate
    // window, so pin/efficiency are grossly mismeasured (a 50 W transformer read
    // as ~150 W, 24.7% "efficiency"). The bound is the physical quantity that
    // caused the bug — the AC-resistance dynamics of a switching magnetic settle
    // within a handful of switching periods (<<100 µs); a stage slower than that
    // is a fit artifact, not physics. The operating-band resistance is unchanged.
    // See abt #71.
    const double maxSettlingTime = 1e-4;  // 100 µs: >> a few switching periods
    for (int k = 0; k < 10; k += 2) {
        const double R = x[k];
        const double L = x[k + 1];
        if (R > 0.0 && L / R > maxSettlingTime) {
            return 1e30;
        }
    }

    return (R0 + parallel(L1, R1 + parallel(L2, R2 + parallel(L3, R3 + parallel(L4, R4 + parallel(L5, R5)))))).real();
}

void CircuitSimulatorExporter::ladder_func(double *p, double *x, int m, int n, void *data) {
    double* dcResistanceAndFrequencies = static_cast <double*> (data);
    double dcResistance = dcResistanceAndFrequencies[0];

    for(int i=0; i<n; ++i)
        x[i]=ladder_model(p, dcResistanceAndFrequencies[i + 1], dcResistance);
}


double CircuitSimulatorExporter::core_ladder_model(double x[], double frequency, double dcResistance) {
    double w = 2 * std::numbers::pi * frequency;
    auto R0 = std::complex<double>(dcResistance);
    auto R1 = std::complex<double>(x[0], 0);
    auto L1 = std::complex<double>(0, w * x[1]);
    auto R2 = std::complex<double>(x[2], 0);
    auto L2 = std::complex<double>(0, w * x[3]);
    auto R3 = std::complex<double>(x[4], 0);
    auto L3 = std::complex<double>(0, w * x[5]);

    // FIXED: Return large penalty instead of 0 to guide optimizer correctly
    for(int i=0; i<6; ++i) {
        if (x[i] < 0) {
            return 1e30;  // Large penalty value instead of 0
        }
    }

    return (R0 + parallel(L1, R1 + parallel(L2, R2 + parallel(L3, R3)))).real();
}

void CircuitSimulatorExporter::core_ladder_func(double *p, double *x, int m, int n, void *data) {
    double* dcResistanceAndFrequencies = static_cast <double*> (data);
    double dcResistance = dcResistanceAndFrequencies[0];

    for(int i=0; i<n; ++i)
        x[i]=core_ladder_model(p, dcResistanceAndFrequencies[i + 1], dcResistance);
}

// Rosano model: R + (L||R) + (L||(R + 1/sC)) — three series stages
// Stage 1 (R):   R1 — pure DC loss floor
// Stage 2 (RL):  L1 || R2 — at DC L shorts→0, at high freq→R2
// Stage 3 (RLC): L2 || (R3 + 1/(sC1)) — adds C for fitting flexibility
// x[0]=R1, x[1]=R2, x[2]=L1, x[3]=R3, x[4]=L2, x[5]=C1
double CircuitSimulatorExporter::core_rosano_model(double x[], double frequency) {
    for (int i = 0; i < 6; ++i) {
        if (x[i] <= 0) {
            return 1e30;
        }
    }

    double w = 2 * std::numbers::pi * frequency;

    // Stage 1: pure R
    auto Z1 = std::complex<double>(x[0], 0);

    // Stage 2: L1 || R2
    auto L1 = std::complex<double>(0, w * x[2]);
    auto R2 = std::complex<double>(x[1], 0);
    auto Z2 = parallel(L1, R2);

    // Stage 3: L2 || (R3 + 1/(jωC1))
    auto L2 = std::complex<double>(0, w * x[4]);
    auto R3plusC = std::complex<double>(x[3], -1.0 / (w * x[5]));
    auto Z3 = parallel(L2, R3plusC);

    // Three stages in series
    return (Z1 + Z2 + Z3).real();
}

// Fitting function in log space for better convergence across wide dynamic range
void CircuitSimulatorExporter::core_rosano_func(double *p, double *x, int m, int n, void *data) {
    double* frequencies = static_cast<double*>(data);

    for (int i = 0; i < n; ++i) {
        double val = core_rosano_model(p, frequencies[i]);
        x[i] = (val > 0) ? std::log(val) : -50.0;
    }
}

// Rosano winding model: Z = Rdc + sum_k(Rk || jwLk)
// Each stage is a resistor in parallel with an inductor, chained in series.
// At DC: all L short → Z = Rdc. At high freq: all L open → Z = Rdc + ΣRk.
// coefficients: [R1, L1, R2, L2, ..., R6, L6]  (12 unknowns, 6 stages)
// Rdc is passed separately (pinned to measured DC resistance, not a fitting parameter).
double CircuitSimulatorExporter::winding_rosano_model(double x[], double frequency, double dcResistance) {
    const int nStages = 6;
    for (int i = 0; i < 2 * nStages; ++i)
        if (x[i] <= 0) return 1e30;

    double w = 2 * std::numbers::pi * frequency;
    auto Z = std::complex<double>(dcResistance, 0);
    for (int k = 0; k < nStages; ++k) {
        auto Rk = std::complex<double>(x[2 * k], 0);
        auto Lk = std::complex<double>(0, w * x[2 * k + 1]);
        Z += parallel(Rk, Lk);
    }
    return Z.real();
}

// Fitting function in log space for uniform relative error across wide dynamic range
void CircuitSimulatorExporter::winding_rosano_func(double *p, double *x, int m, int n, void *data) {
    double* dcResistanceAndFrequencies = static_cast<double*>(data);
    double dcResistance = dcResistanceAndFrequencies[0];
    for (int i = 0; i < n; ++i) {
        double val = winding_rosano_model(p, dcResistanceAndFrequencies[i + 1], dcResistance);
        x[i] = (val > 0) ? std::log(val) : -50.0;
    }
}

// Rosano RLC winding model: Z = Rdc + sum_{k=0}^{nRL-1}(Rk || jwLk) + (Rn || (jwLn + 1/(jwCn)))
// Extended model with one RLC stage for better high-frequency fitting.
// nStages: total number of stages
// nRLCStages: number of RLC stages (1 = one RLC stage at the end, 0 = all RL stages)
// coefficients layout:
//   - If nRLCStages=0: [R1, L1, R2, L2, ..., Rn, Ln] (2*nStages params)
//   - If nRLCStages=1: [R1, L1, ..., R(n-1), L(n-1), Rn, Ln, Cn] (2*nStages+1 params)
// Rdc is passed separately.
double CircuitSimulatorExporter::winding_rosano_rlc_model(double x[], int nStages, int nRLCStages, double frequency, double dcResistance) {
    if (nStages < 1) return dcResistance;
    if (nRLCStages < 0) nRLCStages = 0;
    if (nRLCStages > 1) nRLCStages = 1;

    int nRLStages = nStages - nRLCStages;
    int nParams = 2 * nRLStages + 3 * nRLCStages;

    for (int i = 0; i < nParams; ++i)
        if (x[i] <= 0) return 1e30;

    double w = 2 * std::numbers::pi * frequency;
    auto Z = std::complex<double>(dcResistance, 0);

    for (int k = 0; k < nRLStages; ++k) {
        auto Rk = std::complex<double>(x[2 * k], 0);
        auto Lk = std::complex<double>(0, w * x[2 * k + 1]);
        Z += parallel(Rk, Lk);
    }

    if (nRLCStages > 0) {
        int rlcBase = 2 * nRLStages;
        double R_n = x[rlcBase];
        double L_n = x[rlcBase + 1];
        double C_n = x[rlcBase + 2];

        auto R_branch = std::complex<double>(R_n, 0);
        auto LC_series = std::complex<double>(0, w * L_n - 1.0 / (w * C_n));
        Z += parallel(R_branch, LC_series);
    }

    return Z.real();
}

typedef struct {
    double* dcResistanceAndFrequencies;
    int nStages;
    int nRLCStages;
} WindingRosanoRLCFitData;

void CircuitSimulatorExporter::winding_rosano_rlc_func(double *p, double *x, int m, int n, void *data) {
    auto* fitData = static_cast<WindingRosanoRLCFitData*>(data);
    double dcResistance = fitData->dcResistanceAndFrequencies[0];

    for (int i = 0; i < n; ++i) {
        double val = CircuitSimulatorExporter::winding_rosano_rlc_model(p, fitData->nStages, fitData->nRLCStages,
                                               fitData->dcResistanceAndFrequencies[i + 1], dcResistance);
        x[i] = (val > 0) ? std::log(val) : -50.0;
    }
}

std::vector<std::vector<double>> calculate_ac_resistance_coefficients_per_winding_rosano(Magnetic magnetic, double temperature, bool useRLCStages) {
    const size_t nStages = 6;
    const size_t numberUnknownsRL = 2 * nStages;       // [R1,L1, R2,L2, ..., R6,L6]
    const size_t numberUnknownsRLC = 2 * nStages + 1;  // [R1,L1, ..., R5,L5, R6,L6,C6] - one RLC stage
    const size_t numberElements = 40;
    const size_t numberElementsPlusOne = numberElements + 1;
    double startingFrequency = 0.1;
    double endingFrequency = 10000000;
    auto coil = magnetic.get_coil();

    std::vector<std::vector<double>> acResistanceCoefficientsPerWinding;
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        Curve2D windingAcResistanceData = Sweeper().sweep_winding_resistance_over_frequency(
            magnetic, startingFrequency, endingFrequency, numberElements, windingIndex, temperature);
        auto frequenciesVector = windingAcResistanceData.get_x_points();
        auto valuePoints = windingAcResistanceData.get_y_points();

        double acResistances[numberElements];
        size_t numPoints = std::min(valuePoints.size(), numberElements);
        for (size_t index = 0; index < numPoints; ++index)
            acResistances[index] = valuePoints[index];

        // Pack [Rdc, f0, f1, ...] into data array for fitting function
        double logAcResistances[numberElements];
        double dcResistanceAndFrequencies[numberElementsPlusOne];
        double Rdc = acResistances[0];
        dcResistanceAndFrequencies[0] = Rdc;
        size_t numFreqs = std::min(frequenciesVector.size(), numberElements);
        for (size_t index = 0; index < numFreqs; ++index) {
            dcResistanceAndFrequencies[index + 1] = frequenciesVector[index];
            logAcResistances[index] = std::log(acResistances[index]);
        }

        double Rmax = *std::max_element(valuePoints.begin(), valuePoints.begin() + numPoints);
        double fStart = frequenciesVector.front();
        double fEnd = frequenciesVector.back();

        double bestError = DBL_MAX;
        std::vector<double> bestCoeffs;
        bool bestUseRLC = false;

        // Grid search over rScale and fScale for physics-informed initialization
        double rScales[] = {0.5, 1.0, 2.0, 5.0};
        double fScales[] = {0.3, 1.0, 3.0};

        // First fit RL-only model (standard Rosano)
        for (double rScale : rScales) {
            if (bestError < 0.02) break;
            for (double fScale : fScales) {
                if (bestError < 0.02) break;

                double coefficients[numberUnknownsRL];
                double Rrange = std::max(Rmax - Rdc, Rdc * 0.1) * rScale;

                // Distribute 6 transition frequencies geometrically across [fStart, fEnd]
                for (size_t k = 0; k < nStages; ++k) {
                    double fk = fStart * std::pow(fEnd / fStart, (k + 0.5) / static_cast<double>(nStages)) * fScale;
                    double wk = 2 * std::numbers::pi * fk;
                    double Rk = Rrange / static_cast<double>(nStages);
                    coefficients[2 * k]     = std::max(Rk, 1e-9);
                    coefficients[2 * k + 1] = std::max(Rk / wk, 1e-15);
                }

                double opts[5] = {1e-3, 1e-25, 1e-25, 1e-25, 1e-19};
                double info[10];
                eigen_levmar_dif(CircuitSimulatorExporter::winding_rosano_func, coefficients,
                                 logAcResistances, numberUnknownsRL, numberElements,
                                 10000, opts, info, nullptr, nullptr,
                                 static_cast<void*>(&dcResistanceAndFrequencies));

                bool allPositive = true;
                for (size_t i = 0; i < numberUnknownsRL; ++i)
                    if (coefficients[i] <= 0) { allPositive = false; break; }
                if (!allPositive) continue;

                double errorAverage = 0;
                for (size_t index = 0; index < numFreqs; ++index) {
                    double modeled = CircuitSimulatorExporter::winding_rosano_model(
                        coefficients, frequenciesVector[index], Rdc);
                    double error = fabs(valuePoints[index] - modeled) / valuePoints[index];
                    errorAverage += error;
                }
                errorAverage /= numFreqs;

                if (errorAverage < bestError) {
                    bestError = errorAverage;
                    bestCoeffs = std::vector<double>(coefficients, coefficients + numberUnknownsRL);
                    bestUseRLC = false;
                }
            }
        }

        // Optionally fit RLC model for better high-frequency accuracy
        if (useRLCStages && bestError > 0.01) {  // Only try RLC if RL fit isn't excellent
            WindingRosanoRLCFitData rlcFitData;
            rlcFitData.dcResistanceAndFrequencies = dcResistanceAndFrequencies;
            rlcFitData.nStages = static_cast<int>(nStages);
            rlcFitData.nRLCStages = 1;

            for (double rScale : rScales) {
                if (bestError < 0.01) break;
                for (double fScale : fScales) {
                    if (bestError < 0.01) break;

                    double coefficients[numberUnknownsRLC];
                    double Rrange = std::max(Rmax - Rdc, Rdc * 0.1) * rScale;

                    // 5 RL stages + 1 RLC stage
                    for (size_t k = 0; k < nStages - 1; ++k) {
                        double fk = fStart * std::pow(fEnd / fStart, (k + 0.5) / static_cast<double>(nStages)) * fScale;
                        double wk = 2 * std::numbers::pi * fk;
                        double Rk = Rrange / static_cast<double>(nStages);
                        coefficients[2 * k]     = std::max(Rk, 1e-9);
                        coefficients[2 * k + 1] = std::max(Rk / wk, 1e-15);
                    }
                    // Last stage (RLC)
                    double fk_last = fStart * std::pow(fEnd / fStart, (nStages - 0.5) / static_cast<double>(nStages)) * fScale;
                    double wk_last = 2 * std::numbers::pi * fk_last;
                    double Rk_last = Rrange / static_cast<double>(nStages);
                    int rlcBase = 2 * (nStages - 1);
                    coefficients[rlcBase]     = std::max(Rk_last, 1e-9);
                    coefficients[rlcBase + 1] = std::max(Rk_last / wk_last, 1e-15);
                    coefficients[rlcBase + 2] = std::max(1.0 / (wk_last * Rk_last), 1e-18);

                    double opts[5] = {1e-3, 1e-25, 1e-25, 1e-25, 1e-19};
                    double info[10];
                    eigen_levmar_dif(CircuitSimulatorExporter::winding_rosano_rlc_func, coefficients,
                                     logAcResistances, numberUnknownsRLC, numberElements,
                                     10000, opts, info, nullptr, nullptr,
                                     static_cast<void*>(&rlcFitData));

                    bool allPositive = true;
                    for (size_t i = 0; i < numberUnknownsRLC; ++i)
                        if (coefficients[i] <= 0) { allPositive = false; break; }
                    if (!allPositive) continue;

                    double errorAverage = 0;
                    for (size_t index = 0; index < numFreqs; ++index) {
                        double modeled = CircuitSimulatorExporter::winding_rosano_rlc_model(
                            coefficients, static_cast<int>(nStages), 1, frequenciesVector[index], Rdc);
                        double error = fabs(valuePoints[index] - modeled) / valuePoints[index];
                        errorAverage += error;
                    }
                    errorAverage /= numFreqs;

                    // Only use RLC if it improves error by at least 20%
                    if (errorAverage < bestError * 0.8) {
                        bestError = errorAverage;
                        bestCoeffs = std::vector<double>(coefficients, coefficients + numberUnknownsRLC);
                        bestUseRLC = true;
                    }
                }
            }
        }

        // Append marker flag indicating which model was used. If the fit failed
        // entirely, keep the vector EMPTY so exporters can detect it and fall
        // back to Rdc-only (a lone marker used to defeat their empty() checks,
        // producing a netlist with no path between P+ and the magnetizing node).
        if (!bestCoeffs.empty()) {
            bestCoeffs.push_back(bestUseRLC ? 1.0 : 0.0);
        }

        acResistanceCoefficientsPerWinding.push_back(bestCoeffs);
    }
    return acResistanceCoefficientsPerWinding;
}

std::vector<std::vector<double>> calculate_ac_resistance_coefficients_per_winding_ladder(Magnetic magnetic, double temperature) {
    const size_t numberUnknowns = 10;

    const size_t numberElements = 40;
    const size_t numberElementsPlusOne = numberElements + 1;
    size_t loopIterations = 15;
    double startingFrequency = 0.1;
    double endingFrequency = 10000000;
    auto coil = magnetic.get_coil();

    std::vector<std::vector<double>> acResistanceCoefficientsPerWinding;
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        Curve2D windingAcResistanceData = Sweeper().sweep_winding_resistance_over_frequency(magnetic, startingFrequency, endingFrequency, numberElements, windingIndex, temperature);
        auto frequenciesVector = windingAcResistanceData.get_x_points();

        auto valuePoints = windingAcResistanceData.get_y_points();
        double acResistances[numberElements];
        size_t numPoints = std::min(valuePoints.size(), numberElements);
        for (size_t index = 0; index < numPoints; ++index) {
            acResistances[index] = valuePoints[index];
        }

        double bestError = DBL_MAX;
        double initialState = 10;
        std::vector<double> bestAcResistanceCoefficientsThisWinding;
        for (size_t loopIndex = 0; loopIndex < loopIterations; ++loopIndex) {
            double coefficients[numberUnknowns];
            for (size_t index = 0; index < numberUnknowns; ++index) {
                coefficients[index] = initialState;
            }
            double opts[5], info[10];

            double lmInitMu = 1e-03;
            double lmStopThresh = 1e-25;
            double lmDiffDelta = 1e-19;

            opts[0]= lmInitMu;
            opts[1]= lmStopThresh;
            opts[2]= lmStopThresh;
            opts[3]= lmStopThresh;
            opts[4]= lmDiffDelta;

            double dcResistanceAndfrequencies[numberElementsPlusOne];
            dcResistanceAndfrequencies[0] = acResistances[0];
            coefficients[1] = 1e-9;
            size_t numFreqs = std::min(frequenciesVector.size(), numberElements);
            for (size_t index = 0; index < numFreqs; ++index) {
                dcResistanceAndfrequencies[index + 1] = frequenciesVector[index];
            }

            eigen_levmar_dif(CircuitSimulatorExporter::ladder_func, coefficients, acResistances, numberUnknowns, numberElements, 10000, opts, info, NULL, NULL, static_cast<void*>(&dcResistanceAndfrequencies));

            double errorAverage = 0;
            for (size_t index = 0; index < frequenciesVector.size(); ++index) {
                double modeledAcResistance = CircuitSimulatorExporter::ladder_model(coefficients, frequenciesVector[index], acResistances[0]);
                double error = fabs(valuePoints[index] - modeledAcResistance) / valuePoints[index];
                errorAverage += error;
            }

            errorAverage /= frequenciesVector.size();
            initialState /= 10;

            if (errorAverage < bestError) {
                bestError = errorAverage;
                bestAcResistanceCoefficientsThisWinding.clear();
                for (auto coefficient : coefficients) {
                    bestAcResistanceCoefficientsThisWinding.push_back(coefficient);
                }
            }
        }
        acResistanceCoefficientsPerWinding.push_back(bestAcResistanceCoefficientsThisWinding);

    }
    return acResistanceCoefficientsPerWinding;
}

std::vector<double> CircuitSimulatorExporter::calculate_core_resistance_coefficients(Magnetic magnetic, double temperature, CoreLossTopology topology) {
    const size_t numberUnknowns = 6;

    const size_t numberElements = 20;
    const size_t numberElementsPlusOne = 21;
    size_t loopIterations = 15;
    double startingFrequency = 1000;
    double endingFrequency = 300000;
    auto coil = magnetic.get_coil();

    std::vector<double> coreResistanceCoefficients;


    Curve2D coreResistanceData = Sweeper().sweep_core_resistance_over_frequency(magnetic, startingFrequency, endingFrequency, numberElements, temperature);
    auto frequenciesVector = coreResistanceData.get_x_points();

    auto valuePoints = coreResistanceData.get_y_points();
    double coreResistances[numberElements];
    size_t numPoints = std::min(valuePoints.size(), numberElements);
    for (size_t index = 0; index < numPoints; ++index) {
        coreResistances[index] = valuePoints[index];
    }

    double bestError = DBL_MAX;
    double initialState = 10;
    std::vector<double> bestCoreResistanceCoefficients;

    if (topology == CoreLossTopology::ROSANO) {
        // Rosano model: Z = R1 + (L1||R2) + (L2||(R3+1/sC1))
        // coefficients: [R1, R2, L1, R3, L2, C1]
        // Fit in log space for uniform relative error across wide dynamic range
        double frequencies[numberElements];
        double logCoreResistances[numberElements];
        size_t numFreqs = std::min(frequenciesVector.size(), numberElements);
        for (size_t index = 0; index < numFreqs; ++index) {
            frequencies[index] = frequenciesVector[index];
            logCoreResistances[index] = std::log(coreResistances[index]);
        }

        double Rmin = coreResistances[0];
        double Rmax = coreResistances[numPoints - 1];
        double fMid = std::sqrt(frequenciesVector.front() * frequenciesVector.back());
        double wMid = 2 * std::numbers::pi * fMid;

        // For series model Z = R1 + (L1||R2) + (L2||(R3+1/sC1)):
        // At DC: Z = R1, so R1 ≈ Rmin
        // At high freq: Z ≈ R1 + R2 + R3, so R2+R3 ≈ Rmax
        double rScales[] = {0.5, 1.0, 2.0};
        double freqScales[] = {0.3, 1.0, 3.0};
        double splitRatios[] = {0.5};

        for (double rScale : rScales) {
            if (bestError < 0.05) break;  // Early exit if already good (<5%)
            for (double freqScale : freqScales) {
                if (bestError < 0.05) break;
                for (double splitRatio : splitRatios) {
                double coefficients[numberUnknowns];
                double wTransition = wMid * freqScale;
                double Rrange = (Rmax - Rmin) * rScale;
                if (Rrange <= 0) Rrange = Rmax * rScale;

                // R1: DC floor
                coefficients[0] = Rmin * rScale;
                if (coefficients[0] <= 0) coefficients[0] = 1e-9;
                // R2: RL stage high-freq contribution
                coefficients[1] = Rrange * splitRatio;
                if (coefficients[1] <= 0) coefficients[1] = 1e-6;
                // L1: transition frequency for RL stage
                coefficients[2] = coefficients[1] / wTransition;
                if (coefficients[2] <= 0) coefficients[2] = 1e-9;
                // R3: RLC stage contribution
                coefficients[3] = Rrange * (1.0 - splitRatio);
                if (coefficients[3] <= 0) coefficients[3] = 1e-6;
                // L2: transition frequency for RLC stage
                coefficients[4] = coefficients[3] / (wTransition * 2);
                if (coefficients[4] <= 0) coefficients[4] = 1e-9;
                // C1: resonance around transition frequency
                coefficients[5] = 1.0 / (wTransition * wTransition * coefficients[4]);
                if (coefficients[5] <= 0) coefficients[5] = 1e-9;

                double opts[5], info[10];
                opts[0] = 1e-03;
                opts[1] = 1e-25;
                opts[2] = 1e-25;
                opts[3] = 1e-25;
                opts[4] = 1e-19;

                eigen_levmar_dif(CircuitSimulatorExporter::core_rosano_func, coefficients, logCoreResistances, numberUnknowns, numberElements, 10000, opts, info, NULL, NULL, static_cast<void*>(&frequencies));

                // Check for negative values (invalid physically)
                bool allPositive = true;
                for (size_t i = 0; i < numberUnknowns; ++i) {
                    if (coefficients[i] <= 0) { allPositive = false; break; }
                }
                if (!allPositive) continue;

                double errorAverage = 0;
                for (size_t index = 0; index < frequenciesVector.size(); ++index) {
                    double modeledResistance = CircuitSimulatorExporter::core_rosano_model(coefficients, frequenciesVector[index]);
                    double error = fabs(valuePoints[index] - modeledResistance) / valuePoints[index];
                    errorAverage += error;
                }
                errorAverage /= frequenciesVector.size();

                if (errorAverage < bestError) {
                    bestError = errorAverage;
                    bestCoreResistanceCoefficients.clear();
                    for (auto coefficient : coefficients) {
                        bestCoreResistanceCoefficients.push_back(coefficient);
                    }
                }
            }
            }
        }

        // Post-fit validation: ensure impedance is well-behaved
        if (!bestCoreResistanceCoefficients.empty()) {
            // Check 1: monotonicity — Z should not drop by more than 50%
            double prevZ = 0;
            bool monotonic = true;
            for (size_t i = 0; i < frequenciesVector.size(); ++i) {
                double z = core_rosano_model(bestCoreResistanceCoefficients.data(), frequenciesVector[i]);
                if (z <= 0 || (i > 0 && z < prevZ * 0.5)) {
                    monotonic = false;
                    break;
                }
                prevZ = z;
            }
            if (!monotonic) {
                // Neutralize C1 to remove resonance
                double fMin = frequenciesVector.front();
                double wMin = 2 * std::numbers::pi * fMin;
                double R3 = bestCoreResistanceCoefficients[3];
                bestCoreResistanceCoefficients[5] = 1.0 / (wMin * R3 * 1000.0);
            }

            // Check 2: impedance at mid-frequency must match the Steinmetz data within 10x
            size_t midIdx = frequenciesVector.size() / 2;
            double zModel = core_rosano_model(bestCoreResistanceCoefficients.data(), frequenciesVector[midIdx]);
            double zData = valuePoints[midIdx];
            if (zModel < zData * 0.1 || zModel > zData * 10.0) {
                // Bad fit — discard and return empty (callers handle empty = no core loss network)
                bestCoreResistanceCoefficients.clear();
            }
        }
    } else {
        // Ridley model: R0 + L1 || (R1 + L2 || (R2 + L3 || R3))
        for (size_t loopIndex = 0; loopIndex < loopIterations; ++loopIndex) {
            double coefficients[numberUnknowns];
            for (size_t index = 0; index < numberUnknowns; ++index) {
                coefficients[index] = initialState;
            }
            double opts[5], info[10];

            opts[0] = 1e-03;
            opts[1] = 1e-25;
            opts[2] = 1e-25;
            opts[3] = 1e-25;
            opts[4] = 1e-19;

            double dcResistanceAndfrequencies[numberElementsPlusOne];
            dcResistanceAndfrequencies[0] = coreResistances[0];
            for (size_t index = 0; index < frequenciesVector.size(); ++index) {
                dcResistanceAndfrequencies[index + 1] = frequenciesVector[index];
            }

            eigen_levmar_dif(CircuitSimulatorExporter::core_ladder_func, coefficients, coreResistances, numberUnknowns, numberElements, 10000, opts, info, NULL, NULL, static_cast<void*>(&dcResistanceAndfrequencies));

            double errorAverage = 0;
            for (size_t index = 0; index < frequenciesVector.size(); ++index) {
                double modeledAcResistance = CircuitSimulatorExporter::core_ladder_model(coefficients, frequenciesVector[index], coreResistances[0]);
                double error = fabs(valuePoints[index] - modeledAcResistance) / valuePoints[index];
                errorAverage += error;
            }

            errorAverage /= frequenciesVector.size();
            initialState /= 10;

            if (errorAverage < bestError) {
                bestError = errorAverage;
                bestCoreResistanceCoefficients.clear();
                for (auto coefficient : coefficients) {
                    bestCoreResistanceCoefficients.push_back(coefficient);
                }
            }
        }
        // The Ridley fit pins the series DC resistance to the measured
        // coreResistances[0], which the emitters need to reproduce the fitted
        // impedance. Append it as a trailing odd element: layout becomes
        // [R1, L1, R2, L2, R3, L3, R0] (odd size => last element is R0).
        if (!bestCoreResistanceCoefficients.empty()) {
            bestCoreResistanceCoefficients.push_back(coreResistances[0]);
        }
    }
    return bestCoreResistanceCoefficients;
}


std::vector<std::vector<double>> calculate_ac_resistance_coefficients_per_winding_analytical(Magnetic magnetic, double temperature) {
    const size_t numberUnknowns = 3;  // analytical_model only uses x[0], x[1], x[2]
    const size_t numberElements = 20;

    double startingFrequency = 0.1;
    double endingFrequency = 1000000;
    auto coil = magnetic.get_coil();

    std::vector<std::vector<double>> acResistanceCoefficientsPerWinding;
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        Curve2D windingAcResistanceData = Sweeper().sweep_winding_resistance_over_frequency(magnetic, startingFrequency, endingFrequency, numberElements, windingIndex, temperature);
        auto frequenciesVector = windingAcResistanceData.get_x_points();

        auto points = windingAcResistanceData.get_y_points();
        double acResistances[numberElements];
        size_t numPoints = std::min(points.size(), numberElements);
        for (size_t index = 0; index < numPoints; ++index) {
            acResistances[index] = points[index];
        }

        double coefficients[numberUnknowns];
        for (size_t index = 0; index < numberUnknowns; ++index) {
            coefficients[index] = 1;
        }
        double opts[5], info[10];

        double lmInitMu = 1e-03;
        double lmStopThresh = 1e-25;
        double lmDiffDelta = 1e-19;

        opts[0]= lmInitMu;
        opts[1]= lmStopThresh;
        opts[2]= lmStopThresh;
        opts[3]= lmStopThresh;
        opts[4]= lmDiffDelta;

        double frequencies[numberElements];
        size_t numFreqs = std::min(frequenciesVector.size(), numberElements);
        for (size_t index = 0; index < numFreqs; ++index) {
            frequencies[index] = frequenciesVector[index];
        }

        eigen_levmar_dif(CircuitSimulatorExporter::analytical_func, coefficients, acResistances, numberUnknowns, numberElements, 10000, opts, info, NULL, NULL, static_cast<void*>(&frequencies));

        acResistanceCoefficientsPerWinding.push_back(std::vector<double>());
        for (auto coefficient : coefficients) {
            acResistanceCoefficientsPerWinding.back().push_back(coefficient);
        }

    }
    return acResistanceCoefficientsPerWinding;
}

// =============================================================================
// FRACPOLE RL SYNTHESIS (Kundert DD-profile)
//
// Converts fitted (h, alpha) into concrete R,L values for a series chain of
// parallel RL stages that approximate Z_skin(f) = h * (2*pi*f)^alpha.
// =============================================================================
static std::vector<double> compute_fracpole_rl_coefficients(
    double h, double alpha, double f0, double f1, double lumpsPerDecade = 1.5)
{
    double decades = log10(f1 / f0);
    int N = std::max(2, static_cast<int>(std::round(lumpsPerDecade * decades)));

    double omega0 = 2 * std::numbers::pi * f0;
    double omega1 = 2 * std::numbers::pi * f1;
    double sigma = pow(omega1 / omega0, 1.0 / (2.0 * N));

    // DD profile: poles at omega0 * sigma^(2k), zeros at omega0 * sigma^(2k+1)
    std::vector<double> poles(N + 1), zeros(N);
    for (int k = 0; k <= N; ++k)
        poles[k] = omega0 * pow(sigma, 2.0 * k);
    for (int k = 0; k < N; ++k)
        zeros[k] = omega0 * pow(sigma, 2.0 * k + 1);

    // G(s) = Z_skin(s)/s = h * s^(alpha-1)
    // Rational approximation: G(s) = K * prod(s+zeros) / prod(s+poles)
    double beta = alpha - 1.0;
    double omega_ref = sqrt(omega0 * omega1);

    double H_num_mag = 1.0, H_den_mag = 1.0;
    for (int j = 0; j < N; ++j)
        H_num_mag *= sqrt(omega_ref * omega_ref + zeros[j] * zeros[j]);
    for (int k = 0; k <= N; ++k)
        H_den_mag *= sqrt(omega_ref * omega_ref + poles[k] * poles[k]);

    double target_mag = h * pow(omega_ref, beta);
    double K = target_mag * H_den_mag / H_num_mag;

    // Partial fraction expansion: G(s) = sum_k A_k / (s + poles[k])
    // Z_skin(s) = s * G(s) = sum_k A_k * s / (s + poles[k])
    // Each term = parallel RL stage: R_k = A_k, L_k = A_k / poles[k]
    std::vector<double> coefficients;
    for (int k = 0; k <= N; ++k) {
        double residue_num = K;
        for (int j = 0; j < N; ++j) {
            residue_num *= (zeros[j] - poles[k]);
        }
        double residue_den = 1.0;
        for (int m = 0; m <= N; ++m) {
            if (m != k)
                residue_den *= (poles[m] - poles[k]);
        }
        double A_k = residue_num / residue_den;

        if (A_k > 1e-30) {
            double R_k = A_k;
            double L_k = A_k / poles[k];
            coefficients.push_back(R_k);
            coefficients.push_back(L_k);
        }
    }

    return coefficients;
}

// Compute fracpole coefficients per winding:
// 1. Fit R(f) = Rdc + h*(2*pi*f)^alpha  (2 unknowns: h, alpha)
// 2. Synthesize RL values analytically from Kundert DD-profile algorithm
std::vector<std::vector<double>> calculate_ac_resistance_coefficients_per_winding_fracpole(
    Magnetic magnetic, double temperature)
{
    const size_t numberUnknowns = 2; // h, alpha
    const size_t numberElements = 20;
    const size_t numberElementsPlusOne = numberElements + 1;
    double startingFrequency = 0.1;
    double endingFrequency = 10000000;
    auto coil = magnetic.get_coil();

    std::vector<std::vector<double>> acResistanceCoefficientsPerWinding;
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        Curve2D windingAcResistanceData = Sweeper().sweep_winding_resistance_over_frequency(
            magnetic, startingFrequency, endingFrequency, numberElements, windingIndex, temperature);
        auto frequenciesVector = windingAcResistanceData.get_x_points();
        auto valuePoints = windingAcResistanceData.get_y_points();

        double acResistances[numberElements];
        size_t numPoints = std::min(valuePoints.size(), numberElements);
        for (size_t index = 0; index < numPoints; ++index) {
            acResistances[index] = valuePoints[index];
        }

        double dcResistanceAndFrequencies[numberElementsPlusOne];
        dcResistanceAndFrequencies[0] = acResistances[0]; // DC resistance
        size_t numFreqs = std::min(frequenciesVector.size(), numberElements);
        for (size_t index = 0; index < numFreqs; ++index) {
            dcResistanceAndFrequencies[index + 1] = frequenciesVector[index];
        }

        double coefficients[2];
        coefficients[0] = 1e-4; // Initial guess for h
        coefficients[1] = 0.5;  // Initial guess for alpha (skin effect)

        double opts[5], info[10];
        opts[0] = 1e-03;
        opts[1] = 1e-25;
        opts[2] = 1e-25;
        opts[3] = 1e-25;
        opts[4] = 1e-19;

        eigen_levmar_dif(CircuitSimulatorExporter::fracpole_skin_func,
                         coefficients, acResistances, numberUnknowns, numberElements,
                         10000, opts, info, NULL, NULL,
                         static_cast<void*>(&dcResistanceAndFrequencies));

        double h_fit = std::max(coefficients[0], 1e-15);
        double alpha_fit = std::clamp(coefficients[1], 0.1, 1.0);

        auto rlCoefficients = compute_fracpole_rl_coefficients(
            h_fit, alpha_fit, startingFrequency, endingFrequency, 1.5);

        acResistanceCoefficientsPerWinding.push_back(rlCoefficients);
    }

    return acResistanceCoefficientsPerWinding;
}


std::vector<std::vector<double>> CircuitSimulatorExporter::calculate_ac_resistance_coefficients_per_winding(Magnetic magnetic, double temperature, CircuitSimulatorExporterCurveFittingModes mode) {
    if (mode == CircuitSimulatorExporterCurveFittingModes::LADDER) {
        return calculate_ac_resistance_coefficients_per_winding_ladder(magnetic, temperature);
    }
    else if (mode == CircuitSimulatorExporterCurveFittingModes::FRACPOLE) {
        return calculate_ac_resistance_coefficients_per_winding_fracpole(magnetic, temperature);
    }
    else if (mode == CircuitSimulatorExporterCurveFittingModes::ROSANO) {
        return calculate_ac_resistance_coefficients_per_winding_rosano(magnetic, temperature, false);
    }
    else if (mode == CircuitSimulatorExporterCurveFittingModes::ROSANO_RLC) {
        // Rosano + one terminal RLC stage. The fitter only commits the
        // RLC topology when it improves the average error by ≥ 20 % vs
        // the pure-RL fit, so this mode is a strict superset: worst case
        // it returns identical coefficients to ROSANO.
        return calculate_ac_resistance_coefficients_per_winding_rosano(magnetic, temperature, true);
    }
    else {
        return calculate_ac_resistance_coefficients_per_winding_analytical(magnetic, temperature);
    }
}

// =============================================================================
// MUTUAL RESISTANCE LADDER FITTING
// Based on Hesterman (2020) "Mutual Resistance" paper.
//
// The mutual resistance R_ij represents cross-coupling losses between windings.
// It can be positive (reduces losses when currents oppose) or negative
// (increases losses). The sign typically changes with frequency.
//
// We fit an R-L ladder network to model the frequency dependence:
//   R_mutual(f) = R0 + parallel(L1, R1 + parallel(L2, R2 + ...))
//
// This is similar to the self-resistance ladder but R0 can be negative.
// =============================================================================

// Ladder model for mutual resistance - allows negative values
static double mutual_resistance_ladder_model(double x[], double frequency, double dcResistance) {
    double w = 2 * std::numbers::pi * frequency;
    auto R0 = std::complex<double>(dcResistance);
    auto R1 = std::complex<double>(x[0], 0);
    auto L1 = std::complex<double>(0, w * x[1]);
    auto R2 = std::complex<double>(x[2], 0);
    auto L2 = std::complex<double>(0, w * x[3]);
    auto R3 = std::complex<double>(x[4], 0);
    auto L3 = std::complex<double>(0, w * x[5]);

    // For mutual resistance, parameters can be negative (except inductances)
    // Only check inductances are positive
    if (x[1] < 0 || x[3] < 0 || x[5] < 0) {
        return 1e30;  // Large penalty
    }

    return (R0 + parallel(L1, R1 + parallel(L2, R2 + parallel(L3, R3)))).real();
}

// Data structure for mutual resistance fitting
struct MutualResistanceFitData {
    double dcMutualResistance;
    std::vector<double> frequencies;
};

static void mutual_resistance_ladder_func(double *p, double *x, int m, int n, void *data) {
    MutualResistanceFitData* fitData = static_cast<MutualResistanceFitData*>(data);
    double dcResistance = fitData->dcMutualResistance;
    
    for(int i = 0; i < n; ++i) {
        x[i] = mutual_resistance_ladder_model(p, fitData->frequencies[i], dcResistance);
    }
}

std::vector<CircuitSimulatorExporter::MutualResistanceCoefficients> 
CircuitSimulatorExporter::calculate_mutual_resistance_coefficients(Magnetic magnetic, double temperature) {
    const size_t numberUnknowns = 6;  // 3 R-L pairs for mutual resistance
    const size_t numberElements = 20;
    const size_t loopIterations = 10;
    double startingFrequency = 100;
    double endingFrequency = 1000000;
    
    auto coil = magnetic.get_coil();
    size_t numWindings = coil.get_functional_description().size();
    
    std::vector<MutualResistanceCoefficients> allMutualCoeffs;
    
    // Only calculate for multi-winding components
    if (numWindings < 2) {
        return allMutualCoeffs;
    }
    
    // Get resistance matrix at multiple frequencies
    auto resistanceMatrices = Sweeper::sweep_resistance_matrix_over_frequency(
        magnetic, startingFrequency, endingFrequency, numberElements, temperature, "log");
    
    if (resistanceMatrices.empty()) {
        return allMutualCoeffs;
    }
    
    // Extract frequencies
    std::vector<double> frequencies;
    for (const auto& matrix : resistanceMatrices) {
        frequencies.push_back(matrix.get_frequency());
    }
    
    // Get winding names
    auto& functionalDescription = coil.get_functional_description();
    std::vector<std::string> windingNames;
    for (const auto& winding : functionalDescription) {
        windingNames.push_back(winding.get_name());
    }
    
    // Calculate coefficients for each winding pair
    for (size_t i = 0; i < numWindings; ++i) {
        for (size_t j = i + 1; j < numWindings; ++j) {
            MutualResistanceCoefficients mutualCoeffs;
            mutualCoeffs.windingIndex1 = i;
            mutualCoeffs.windingIndex2 = j;
            
            // Extract mutual resistance values at each frequency
            std::vector<double> mutualResistances(numberElements);
            for (size_t k = 0; k < numberElements; ++k) {
                auto& magnitude = resistanceMatrices[k].get_magnitude();
                auto it_i = magnitude.find(windingNames[i]);
                if (it_i != magnitude.end()) {
                    auto it_j = it_i->second.find(windingNames[j]);
                    if (it_j != it_i->second.end()) {
                        mutualResistances[k] = it_j->second.get_nominal().value();
                    }
                }
            }
            
            // DC mutual resistance (from lowest frequency)
            mutualCoeffs.dcMutualResistance = mutualResistances[0];
            
            // Prepare data for fitting
            MutualResistanceFitData fitData;
            fitData.dcMutualResistance = mutualCoeffs.dcMutualResistance;
            fitData.frequencies = frequencies;
            
            double bestError = DBL_MAX;
            double initialState = 1.0;
            std::vector<double> bestCoeffs;
            
            for (size_t loopIndex = 0; loopIndex < loopIterations; ++loopIndex) {
                double coefficients[numberUnknowns];
                for (size_t idx = 0; idx < numberUnknowns; ++idx) {
                    coefficients[idx] = initialState;
                }
                // Initialize inductances to small positive values
                coefficients[1] = 1e-7;
                coefficients[3] = 1e-8;
                coefficients[5] = 1e-9;
                
                double opts[5], info[10];
                opts[0] = 1e-03;
                opts[1] = 1e-20;
                opts[2] = 1e-20;
                opts[3] = 1e-20;
                opts[4] = 1e-17;
                
                double targetValues[numberElements];
                for (size_t idx = 0; idx < numberElements; ++idx) {
                    targetValues[idx] = mutualResistances[idx];
                }
                
                eigen_levmar_dif(mutual_resistance_ladder_func, coefficients, targetValues, 
                                 numberUnknowns, numberElements, 5000, opts, info, NULL, NULL,
                                 static_cast<void*>(&fitData));
                
                // Calculate fitting error
                double errorAverage = 0;
                for (size_t idx = 0; idx < numberElements; ++idx) {
                    double modeled = mutual_resistance_ladder_model(coefficients, frequencies[idx], 
                                                                    mutualCoeffs.dcMutualResistance);
                    double target = mutualResistances[idx];
                    // Use relative error for non-zero values, absolute for near-zero
                    if (std::abs(target) > 1e-6) {
                        errorAverage += std::abs(target - modeled) / std::abs(target);
                    } else {
                        errorAverage += std::abs(target - modeled);
                    }
                }
                errorAverage /= numberElements;
                initialState /= 10;
                
                if (errorAverage < bestError) {
                    bestError = errorAverage;
                    bestCoeffs.clear();
                    for (size_t idx = 0; idx < numberUnknowns; ++idx) {
                        bestCoeffs.push_back(coefficients[idx]);
                    }
                }
            }
            
            mutualCoeffs.coefficients = bestCoeffs;
            allMutualCoeffs.push_back(mutualCoeffs);
        }
    }
    
    return allMutualCoeffs;
}


std::shared_ptr<CircuitSimulatorExporterModel> CircuitSimulatorExporterModel::factory(CircuitSimulatorExporterModels programName){
    if (programName == CircuitSimulatorExporterModels::SIMBA) {
        return std::make_shared<CircuitSimulatorExporterSimbaModel>();
    }
    else if (programName == CircuitSimulatorExporterModels::NGSPICE) {
        return std::make_shared<CircuitSimulatorExporterNgspiceModel>();
    }
    else if (programName == CircuitSimulatorExporterModels::LTSPICE) {
        return std::make_shared<CircuitSimulatorExporterLtspiceModel>();
    }
    else if (programName == CircuitSimulatorExporterModels::PLECS) {
        return std::make_shared<CircuitSimulatorExporterPlecsModel>();
    }
    else if (programName == CircuitSimulatorExporterModels::NL5) {
        // ABT #120.1: was missing — NL5 is in the enum with a full exporter class,
        // but the factory threw "Unknown program" for it.
        return std::make_shared<CircuitSimulatorExporterNl5Model>();
    }
    else
        throw ModelNotAvailableException("Unknown Circuit Simulator program, available options are: {SIMBA, NGSPICE, LTSPICE, PLECS, NL5}");
}

std::string CircuitSimulatorExporter::export_magnetic_as_subcircuit(Magnetic magnetic, double frequency, double temperature, std::optional<std::string> outputFilename, std::optional<std::string> filePathOrFile, CircuitSimulatorExporterCurveFittingModes mode) {
    auto result = _model->export_magnetic_as_subcircuit(magnetic, frequency, temperature, filePathOrFile, mode);
    if (outputFilename) {
        std::ofstream o(outputFilename.value());
        o << std::setw(2) << result << std::endl;
    }
    return result;
}

std::string CircuitSimulatorExporter::export_magnetic_as_symbol(Magnetic magnetic, std::optional<std::string> outputFilename, std::optional<std::string> filePathOrFile) {
    auto result = _model->export_magnetic_as_symbol(magnetic, filePathOrFile);
    if (outputFilename) {
        std::ofstream o(outputFilename.value());
        o << std::setw(2) << result << std::endl;
    }
    return result;
}


// =============================================================================
// FRACPOLE SPICE EMITTERS
// =============================================================================

// Emit the fracpole subcircuit + instantiation for winding i in NgSpice/LTSpice format.
// The winding resistance model is:
//   Rdc in series with the skin_effect subcircuit (fracpole + gyrator).
//   P+ ---[Rdc]---[skin_effect]--- Node_R_Lmag
std::string emit_fracpole_winding_spice(
    const FractionalPoleNetwork& net, size_t windingIndex, const std::string& is,
    double dcResistance) {
    std::string s;
    const std::string fp_name = "fracpole_w" + is;
    const std::string sk_name = "skin_effect_w" + is;
    double coef = net.opts.coef;

    // Fracpole subcircuit definition (R/C stages) - pre-computed values for ngspice compatibility
    s += "* Winding " + is + " skin-effect fracpole (alpha=0.5, "
       + std::to_string(net.stages.size()) + " stages)\n";
    s += ".subckt " + fp_name + " p n\n";
    for (size_t k = 0; k < net.stages.size(); ++k) {
        std::string kk = std::to_string(k + 1);
        double R_scaled = net.stages[k].R * coef;
        double C_scaled = net.stages[k].C / coef;
        s += "R" + kk + " p " + kk + " " + to_string(R_scaled, 12) + "\n";
        s += "C" + kk + " " + kk + " n " + to_string(C_scaled, 12) + "\n";
    }
    if (net.Cinf > 0.0) {
        double Cinf_scaled = net.Cinf / coef;
        s += "Cinf p n " + to_string(Cinf_scaled, 12) + "\n";
    }
    if (net.Rinf > 0.0) {
        double Rinf_scaled = net.Rinf * coef;
        s += "Rinf p n " + to_string(Rinf_scaled, 12) + "\n";
    }
    s += ".ends " + fp_name + "\n\n";

    // Skin-effect subcircuit definition (fracpole + gyrator)
    // The gyrator converts impedance fracpole -> admittance fracpole
    // so resistance increases with frequency (skin effect)
    s += ".subckt " + sk_name + " 1 2\n";
    s += "* Gyrator: converts impedance fracpole to admittance fracpole\n";
    s += "Eg1 3 0 1 2 1\n";          // VCVS sense: input port voltage
    s += "Gg1 1 2 3 0 1\n";          // VCCS: drives port with gyrator current
    s += "Rg1 3 0 1e12\n";           // High impedance termination (gmin^-1)
    s += "Xfp 3 0 " + fp_name + "\n";
    s += ".ends " + sk_name + "\n\n";

    // DC resistance element
    s += "Rdc" + is + " P" + is + "+ Node_fp_" + is + " {Rdc_" + is + "_Value}\n";
    // Skin effect instance
    s += "Xskin" + is + " Node_fp_" + is + " Node_R_Lmag_" + is + " " + sk_name + "\n";

    return s;
}

// Emit the fracpole subcircuit for core losses.
// Core loss is represented as a frequency-dependent resistance in parallel
// with the magnetizing inductance.
// Uses pre-computed values (coef already applied) for ngspice compatibility.
//
// IMPORTANT: For core losses, we do NOT include Cinf or Rinf terminal elements.
// - Cinf would short-circuit the magnetizing inductance at high frequencies
// - Rinf would provide a DC path that doesn't exist for core losses
// The R-C stages alone provide the frequency-dependent loss behavior we need.
// NOTE: emit_fracpole_core_spice is currently unused; kept for reference.
#if 0
static std::string emit_fracpole_core_spice(
    const FractionalPoleNetwork& net, size_t numWindings) {
    std::string s;
    const std::string fp_name = "fracpole_core";
    double coef = net.opts.coef;

    s += "* Core loss fracpole (alpha=" + to_string(net.opts.alpha, 3) + ", "
       + std::to_string(net.stages.size()) + " stages)\n";
    s += ".subckt " + fp_name + " p n\n";
    for (size_t k = 0; k < net.stages.size(); ++k) {
        std::string kk = std::to_string(k + 1);
        // Apply coef scaling directly to values for ngspice compatibility
        double R_scaled = net.stages[k].R * coef;
        double C_scaled = net.stages[k].C / coef;
        s += "R" + kk + " p " + kk + " " + to_string(R_scaled, 12) + "\n";
        s += "C" + kk + " " + kk + " n " + to_string(C_scaled, 12) + "\n";
    }
    // NOTE: Cinf and Rinf are intentionally omitted for core losses.
    // Cinf would short-circuit Lmag at high frequencies, causing simulation failure.
    // Rinf would provide an unphysical DC path through the core loss network.
    s += ".ends " + fp_name + "\n\n";

    // Instantiate in parallel with the magnetizing inductance node
    s += "Xcore_loss Node_R_Lmag_1 P1- " + fp_name + "\n";

    return s;
}

#endif
// Emit core loss network using R-L ladder topology for SPICE.
// The R-L ladder provides impedance that INCREASES with frequency:
// - At low frequency: L shorts out R, so total Z is low
// - At high frequency: L is open, so R dominates, Z is high
// This matches the behavior needed for core losses (higher losses at higher frequency).
// Coefficients are [R1, L1, R2, L2, R3, L3] from calculate_core_resistance_coefficients(),
// plus a trailing R0 (measured DC core resistance) when the size is odd.
// Fitted model (core_ladder_model): Z = R0 + L1 || (R1 + L2 || (R2 + L3 || R3))
std::string emit_core_ladder_spice(
    const std::vector<double>& coeffs, size_t numWindings) {

    if (coeffs.size() < 2) {
        return "";  // Not enough coefficients
    }

    double dcResistance = 0;
    size_t stageCount = coeffs.size() / 2;
    if (coeffs.size() % 2 == 1) {
        dcResistance = coeffs.back();
    }

    // Collect the leading run of valid stages; a broken stage terminates the
    // ladder instead of being silently skipped (skipping would change the
    // topology of the remaining stages).
    std::vector<std::pair<double, double>> stages;
    for (size_t k = 0; k < stageCount; ++k) {
        double R = coeffs[k * 2];
        double L = coeffs[k * 2 + 1];
        if (R <= 0 || L <= 0 || R > 1e6 || L > 1) {
            break;
        }
        stages.push_back({R, L});
    }

    if (stages.empty() && dcResistance <= 0) {
        return "";
    }

    std::string s;
    s += "* Core loss R-L ladder (" + std::to_string(stages.size()) + " stages): Z = R0 + L1||(R1 + L2||(R2 + ...))\n";

    // Series R0 from the magnetizing node into the ladder
    std::string currentNode = "Node_R_Lmag_1";
    if (dcResistance > 0) {
        std::string firstLadderNode = stages.empty() ? "P1-" : "Node_Rcore_0";
        s += "Rcore0 " + currentNode + " " + firstLadderNode + " " + to_string(dcResistance, 12) + "\n";
        currentNode = firstLadderNode;
    }

    for (size_t k = 0; k < stages.size(); ++k) {
        auto [R, L] = stages[k];
        std::string kk = std::to_string(k + 1);
        std::string nextNode = (k + 1 == stages.size()) ? "P1-" : ("Node_Rcore_" + std::to_string(k + 1));
        // L_k in parallel with (R_k + rest): L from this node to P1-
        s += "Lcore" + kk + " " + currentNode + " P1- " + to_string(L, 12) + "\n";
        // R_k in the series path towards the next stage (or P1- for the last)
        s += "Rcore" + kk + " " + currentNode + " " + nextNode + " " + to_string(R, 12) + "\n";
        currentNode = nextNode;
    }

    return s;
}

// Emit core loss network using Rosano topology for SPICE.
// Three stages in series from Node_R_Lmag_1 to P1-:
//   Stage 1 (R):   R1 — pure resistance
//   Stage 2 (RL):  L1 || R2 — inductor in parallel with resistor
//   Stage 3 (RLC): L2 || (R3 + C1) — inductor in parallel with resistor-capacitor series
// Coefficients are [R1, R2, L1, R3, L2, C1] from calculate_core_resistance_coefficients(ROSANO).
std::string emit_core_rosano_spice(
    const std::vector<double>& coeffs, size_t numWindings) {

    if (coeffs.size() < 6) {
        return "";
    }

    double R1 = coeffs[0];
    double R2 = coeffs[1];
    double L1 = coeffs[2];
    double R3 = coeffs[3];
    double L2 = coeffs[4];
    double C1 = coeffs[5];

    // The network can only be terminated to P1- (ground) through the RLC
    // stage 3 (L2 || (R3+C1)). Every coefficient must be finite and positive
    // for the R + RL + RLC topology to be well-formed. If the fit is
    // degenerate, DO NOT close the branch with the old `Rcore_gnd .. 1e-6`
    // fallback: with stages 1/2 also skipped that dropped a ~1 uOhm resistor
    // straight across Node_R_Lmag_1..P1- and SHORTED the magnetizing
    // inductance (ABT #120.3). Return "" so the caller simply omits the
    // core-loss branch (Lmag stands alone) rather than shorting it.
    if (!std::isfinite(R1) || !std::isfinite(R2) || !std::isfinite(L1) ||
        !std::isfinite(R3) || !std::isfinite(L2) || !std::isfinite(C1) ||
        R1 <= 0 || R2 <= 0 || L1 <= 0 || R3 <= 0 || L2 <= 0 || C1 <= 0) {
        return "";
    }

    std::string s;
    s += "* Core loss Rosano network (R + RL + RLC series stages)\n";

    // Stage 1 (R): series R from core node
    std::string prevNode = "Node_R_Lmag_1";
    std::string nextNode = "Node_core_s1";
    s += "Rcore_s1 " + prevNode + " " + nextNode + " " + to_string(R1, 12) + "\n";
    prevNode = nextNode;

    // Stage 2 (RL): L1 || R2 between prevNode and nextNode
    nextNode = "Node_core_s2";
    // R2 path
    s += "Rcore_s2 " + prevNode + " " + nextNode + " " + to_string(R2, 12) + "\n";
    // L1 in parallel (same two nodes)
    s += "Lcore_s2 " + prevNode + " " + nextNode + " " + to_string(L1, 12) + "\n";
    prevNode = nextNode;

    // Stage 3 (RLC): L2 || (R3 + C1) between prevNode and ground
    // L2 path: prevNode to ground
    s += "Lcore_s3 " + prevNode + " P1- " + to_string(L2, 12) + "\n";
    // R3 + C1 path (in series): prevNode → R3 → C1 → ground
    s += "Rcore_s3 " + prevNode + " Node_core_s3c " + to_string(R3, 12) + "\n";
    s += "Ccore_s3 Node_core_s3c P1- " + to_string(C1, 12) + "\n";

    return s;
}

// Emit Rosano winding AC-resistance network for SPICE.
// Topology: Z = Rdc + (R1||L1) + (R2||L2) + ... + (R6||L6)  [flat series of parallel R||L cells]
// Rdc is emitted separately. Each cell Rk||Lk sits between two successive nodes.
// coefficients: [R1, L1, R2, L2, ..., R6, L6]
std::string emit_winding_rosano_spice(
    const std::vector<double>& coeffs,
    const std::string& is,
    double dcResistance) {

    if (coeffs.size() < 2) return "";

    // Check if RLC model was used (marker flag appended as last element)
    bool useRLC = false;
    size_t actualCoeffCount = coeffs.size();
    if (actualCoeffCount > 0 && coeffs.back() == 1.0) {
        useRLC = true;
        actualCoeffCount--;
    } else if (actualCoeffCount > 0 && coeffs.back() == 0.0) {
        actualCoeffCount--;
    }

    size_t nStages;
    if (useRLC) {
        nStages = (actualCoeffCount - 1) / 2;  // 5 RL stages + 1 RLC stage
    } else {
        nStages = actualCoeffCount / 2;
    }

    std::string s;
    s += "* Rosano winding AC loss network (" + std::to_string(nStages) + " stages, " +
         (useRLC ? "RLC" : "RL") + " model)\n";

    // Rdc in series (R_0 in Rosano's notation)
    s += "Rdc" + is + " P" + is + "+ Node_Rdc_" + is + " " + to_string(dcResistance, 12) + "\n";

    std::string prevNode = "Node_Rdc_" + is;
    for (size_t k = 0; k < nStages; ++k) {
        double Rk = coeffs[2 * k];
        double Lk = coeffs[2 * k + 1];
        std::string ks = std::to_string(k + 1);
        bool isLastStage = (k == nStages - 1);
        std::string nextNode = isLastStage ? "Node_R_Lmag_" + is : ("Node_Rw_" + is + "_" + ks);

        if (useRLC && isLastStage) {
            // RLC stage: R || (L + C series)
            double Ck = coeffs[actualCoeffCount - 1];
            if (Rk > 0 && Lk > 0 && Ck > 0) {
                s += "Rw" + is + "_" + ks + " " + prevNode + " " + nextNode + " " + to_string(Rk, 12) + "\n";
                s += "Lw" + is + "_" + ks + " " + prevNode + " Node_Lw_" + is + "_" + ks + " " + to_string(Lk, 12) + "\n";
                s += "Cw" + is + "_" + ks + " Node_Lw_" + is + "_" + ks + " " + nextNode + " " + to_string(Ck, 12) + "\n";
            } else {
                s += "Rwshort" + is + "_" + ks + " " + prevNode + " " + nextNode + " 1e-9\n";
            }
        } else {
            // RL stage: R || L
            if (Rk > 0 && Lk > 0) {
                s += "Rw" + is + "_" + ks + " " + prevNode + " " + nextNode + " " + to_string(Rk, 12) + "\n";
                s += "Lw" + is + "_" + ks + " " + prevNode + " " + nextNode + " " + to_string(Lk, 12) + "\n";
            } else {
                s += "Rwshort" + is + "_" + ks + " " + prevNode + " " + nextNode + " 1e-9\n";
            }
        }
        prevNode = nextNode;
    }
    return s;
}



// =============================================================================
// MUTUAL RESISTANCE AUXILIARY WINDING NETWORK
// Based on Hesterman (2020) "Mutual Resistance" paper.
//
// Implements the mutual resistance model using auxiliary windings coupled
// to the physical windings. Each winding pair has an auxiliary winding
// network that models the cross-coupling AC losses.
//
// Circuit topology for each winding pair (i, j):
//   Physical winding i (Lmag_i) is coupled to auxiliary winding LA_ij
//   Physical winding j (Lmag_j) is coupled to auxiliary winding LA_ij
//   Shunt resistor RA_ij on auxiliary winding dissipates mutual losses
//
// The coupling coefficients determine how the mutual resistance varies
// with frequency based on the fitted R-L ladder coefficients.
// =============================================================================

std::string emit_mutual_resistance_network_spice(
    const std::vector<CircuitSimulatorExporter::MutualResistanceCoefficients>& mutualCoeffs,
    double magnetizingInductance,
    size_t numWindings) {
    
    if (mutualCoeffs.empty() || numWindings < 2) {
        return "";  // No mutual resistance to model
    }

    // The auxiliary-winding mutual-resistance model couples an extra inductor LA_ij (k=0.1)
    // into BOTH magnetizing inductors of every winding pair. For a 2-winding transformer the
    // resulting coupled-inductor system {Lmag_1, Lmag_2, LA_12} is small and positive-definite,
    // so ngspice solves it. For n >= 3 windings the pairs share the magnetizing inductors
    // (LA_12 and LA_13 both couple Lmag_1) but carry NO direct LA-LA coupling, so ngspice sees a
    // single coupled system whose K matrix is INCOMPLETE (the missing entries are implicitly 0)
    // and NOT positive definite -> the transient collapses ("Timestep too small"). The mutual
    // (cross-coupling) resistance is a SECOND-ORDER loss term; the dominant per-winding self AC
    // resistance (Rdc + ladder) is emitted separately and is unaffected. Rather than emit an
    // unsimulatable deck, skip this term for n >= 3 and say so loudly. (A PD-safe behavioural
    // current-controlled implementation that keeps the term for n >= 3 is follow-up — abt #50.)
    if (numWindings >= 3) {
        return "\n* ==== MUTUAL RESISTANCE NETWORK: SKIPPED ====\n"
               "* Skipped for " + std::to_string(numWindings) + "-winding transformer: the "
               "auxiliary-winding (Hesterman 2020) coupled-inductor realization yields a "
               "non-positive-definite\n* coupled-L matrix in ngspice for >=3 windings (incomplete "
               "LA-LA K couplings). The dominant per-winding\n* self AC resistance is still "
               "modelled; only the second-order cross-coupling loss term is omitted (abt #50).\n";
    }

    std::string s;
    s += "\n* ==== MUTUAL RESISTANCE NETWORK ====\n";
    s += "* Based on Hesterman (2020) - Auxiliary windings for cross-coupling losses\n";
    s += "* Power: P = R11*I1^2 + R22*I2^2 + 2*R12*I1*I2*cos(theta)\n\n";
    
    for (const auto& mc : mutualCoeffs) {
        size_t i = mc.windingIndex1 + 1;  // 1-indexed
        size_t j = mc.windingIndex2 + 1;
        std::string ij = std::to_string(i) + std::to_string(j);
        
        // Check if we have valid coefficients
        if (mc.coefficients.size() < 2) {
            continue;
        }
        
        // Determine the coupling polarity based on DC mutual resistance sign
        // Positive R12: reduces losses when currents oppose (typical transformer)
        // Negative R12: increases losses
        double sign = (mc.dcMutualResistance >= 0) ? 1.0 : -1.0;
        
        s += "* Mutual resistance R" + ij + " (dc=" + to_string(mc.dcMutualResistance, 6) + " Ohm)\n";
        
        // Create auxiliary winding network with R-L ladder
        // The auxiliary winding models the frequency-dependent mutual resistance
        // Topology: series of (L || R) stages
        size_t numStages = mc.coefficients.size() / 2;
        
        // First, create the auxiliary winding with the same inductance as Lmag
        // This auxiliary winding is coupled to both physical windings
        s += "LA_" + ij + " Node_LA_" + ij + "_in Node_LA_" + ij + "_out " + to_string(magnetizingInductance, 12) + "\n";
        
        // Add shunt resistor network on auxiliary winding to model losses
        std::string lastNode = "Node_LA_" + ij + "_out";
        for (size_t stage = 0; stage < numStages; ++stage) {
            double R = std::abs(mc.coefficients[stage * 2]);
            double L = mc.coefficients[stage * 2 + 1];
            std::string sk = std::to_string(stage);
            
            // Validate values
            if (R <= 0 || R > 1e6 || L <= 0 || L > 1) {
                continue;
            }
            
            std::string nodeR = "Node_Rmut_" + ij + "_" + sk;
            
            // L in parallel with R-chain
            s += "Lmut_" + ij + "_" + sk + " " + lastNode + " " + nodeR + " " + to_string(L, 12) + "\n";
            // R in series to next stage
            if (stage == numStages - 1) {
                // Last stage connects to ground
                s += "Rmut_" + ij + "_" + sk + " " + nodeR + " 0 " + to_string(R, 12) + "\n";
            } else {
                s += "Rmut_" + ij + "_" + sk + " " + nodeR + " Node_Rmut_" + ij + "_" + std::to_string(stage + 1) + " " + to_string(R, 12) + "\n";
            }
            
            lastNode = nodeR;
        }
        
        // Connect auxiliary winding input to ground (forms the load)
        s += "Rmut_" + ij + "_gnd Node_LA_" + ij + "_in 0 1e9\n";  // High-Z to avoid DC path
        
        // Coupling coefficients for auxiliary winding
        // The auxiliary winding is loosely coupled to both physical windings
        // The coupling determines how current in each winding induces loss in the auxiliary
        double kBase = 0.1;  // Base coupling - small for loose coupling
        
        // Scale by sqrt of ratio to balance energy transfer
        double k_i_aux = sign * kBase;
        double k_j_aux = sign * kBase;
        
        // Clamp coupling to valid range (-1, 1)
        k_i_aux = std::max(-0.99, std::min(0.99, k_i_aux));
        k_j_aux = std::max(-0.99, std::min(0.99, k_j_aux));
        
        s += "KA_" + std::to_string(i) + "_" + ij + " Lmag_" + std::to_string(i) + " LA_" + ij + " " + to_string(k_i_aux, 6) + "\n";
        s += "KA_" + std::to_string(j) + "_" + ij + " Lmag_" + std::to_string(j) + " LA_" + ij + " " + to_string(k_j_aux, 6) + "\n";

        s += "\n";
    }

    return s;
}


std::string emit_mutual_resistance_behavioural_spice(
    const std::vector<CircuitSimulatorExporter::MutualResistanceCoefficients>& mutualCoeffs,
    size_t numWindings,
    double frequency) {

    if (mutualCoeffs.empty() || numWindings < 3) {
        return "";  // Behavioural realization is only needed for n>=3 (n==2 is PD already)
    }

    // PD-safe behavioural realization of the cross-coupling (mutual) resistance loss for
    // n>=3 windings (abt #50/#72). The auxiliary-winding (Hesterman 2020) model couples an
    // extra inductor into TWO magnetizing inductors at once; for n>=3 the shared-but-
    // uncoupled auxiliaries leave the coupled-L matrix non-positive-definite and ngspice
    // collapses ("Timestep too small"), so it was skipped. Here each winding i instead
    // carries a series voltage drop sum_{j!=i} R_ij*I_j realized with LINEAR current-
    // controlled voltage sources (H), so NO inductor and NO coupling is added: the Lmag
    // coupled-L matrix is byte-identical to the no-mutual case and stays positive-definite.
    //
    //   * Vsns_<i>      : a 0 V sense source, so i(Vsns_<i>) = the winding-i terminal current
    //   * Hmut_<i>_<j>  : a CCVS = R_ij * i(Vsns_<j>), chained in series into the drop
    // R_ij is the PURELY RESISTIVE mutual resistance at the export frequency (the real part
    // of the fitted ladder impedance, mutual_resistance_real_at). Using a scalar resistance
    // rather than the raw R-L ladder keeps the loss WITHOUT the ladder's reactance, which
    // would otherwise add a spurious series inductance and detune the winding — the reactive
    // coupling already lives in the Lmag K statements. The winding emission routes
    // P<i>+ -> Vsns_<i> -> H chain -> Node_Wtop_<i> -> Rdc. Power: the H drop in winding i
    // absorbs (R_ij*I_j)*I_i and in winding j absorbs (R_ij*I_i)*I_j — the full
    // 2*R_ij*I_i*I_j Hesterman cross term, at the winding terminals.

    // Per receiver winding (0-indexed): list of (driving winding, signed R_ij) drops.
    std::vector<std::vector<std::pair<size_t, double>>> drops(numWindings);

    std::string s;
    s += "\n* ==== MUTUAL RESISTANCE NETWORK (PD-safe behavioural, n>=3) ====\n";
    s += "* Cross-coupling loss sum_{i!=j} R_ij*I_i*I_j via linear current-controlled\n";
    s += "* voltage sources (no inductors, no coupling) so the Lmag matrix stays\n";
    s += "* positive-definite (abt #50/#72). R_ij is the mutual resistance at the export\n";
    s += "* frequency (real part of the fitted ladder; the reactance lives in the K's).\n\n";

    for (const auto& mc : mutualCoeffs) {
        if (mc.coefficients.size() < 2) {
            continue;
        }
        if (mc.windingIndex1 >= numWindings || mc.windingIndex2 >= numWindings ||
            mc.windingIndex1 == mc.windingIndex2) {
            continue;
        }
        size_t i = mc.windingIndex1;   // 0-indexed
        size_t j = mc.windingIndex2;
        // Mutual resistance (signed) at the export frequency, evaluated with the SAME
        // model the coefficients were fitted to (mutual_resistance_ladder_model). This is
        // a pure scalar resistance carrying the DC term and the correct sign — no ladder,
        // so no spurious reactance (that reactive coupling already lives in the K's).
        if (mc.coefficients.size() < 6) {
            continue;  // the fit always emits 6 coefficients; skip a malformed pair
        }
        double rij = mutual_resistance_ladder_model(
            const_cast<double*>(mc.coefficients.data()), frequency, mc.dcMutualResistance);
        if (!std::isfinite(rij) || rij == 0.0) {
            continue;
        }
        s += "* Mutual R" + std::to_string(i + 1) + std::to_string(j + 1) +
             " = " + to_string(rij, 9) + " Ohm @ " + to_string(frequency, 1) + " Hz\n";
        drops[i].emplace_back(j, rij);
        drops[j].emplace_back(i, rij);
    }

    // Per-winding series sense (Vsns) + a chain of CCVS drops (Hmut), routing P<k>+ into
    // the Node_Wtop_<k> node the winding emission consumes. Every winding gets a
    // Node_Wtop_<k> (a straight 0 V wire when it has no mutual term) so the winding side
    // never references an undefined node.
    for (size_t k = 0; k < numWindings; ++k) {
        std::string ks = std::to_string(k + 1);
        s += "Vsns_" + ks + " P" + ks + "+ Node_Wsns_" + ks + " 0\n";
        const auto& kdrops = drops[k];
        if (kdrops.empty()) {
            s += "Vwire_" + ks + " Node_Wsns_" + ks + " Node_Wtop_" + ks + " 0\n";
            continue;
        }
        std::string node = "Node_Wsns_" + ks;
        for (size_t d = 0; d < kdrops.size(); ++d) {
            std::string otherS = std::to_string(kdrops[d].first + 1);
            std::string next = (d + 1 == kdrops.size())
                                   ? ("Node_Wtop_" + ks)
                                   : ("Node_Hmut_" + ks + "_" + std::to_string(d));
            s += "Hmut_" + ks + "_" + otherS + " " + node + " " + next + " Vsns_" + otherS +
                 " " + to_string(kdrops[d].second, 9) + "\n";
            node = next;
        }
    }
    s += "\n";
    return s;
}

// Stray/parasitic capacitance network (positive 3-capacitor / pi-model). Shared verbatim by the
// ngspice and LTspice exporters — the caps are plain `C…` terminal elements with identical syntax
// in both. See the header for the model rationale.
std::string emit_stray_capacitance_spice(const Coil& coil, size_t numWindings) {
    auto& settings = Settings::GetInstance();
    if (!settings.get_circuit_simulator_include_stray_capacitance() || !coil.get_turns_description()) {
        return "";
    }
    auto capAmongWindings = StrayCapacitance().calculate_capacitance(coil).get_capacitance_among_windings();
    if (!capAmongWindings) {
        return "";
    }
    const auto& capMap = capAmongWindings.value();
    const auto& fd = coil.get_functional_description();
    std::string capString = "\n* ==== STRAY CAPACITANCE NETWORK ====\n";
    capString += "* Per-winding self-capacitance (self-resonance) + inter-winding\n";
    capString += "* capacitance (CM coupling), lumped positive caps from StrayCapacitance.\n";
    bool emittedAny = false;
    for (size_t a = 0; a < numWindings; ++a) {
        std::string an = fd[a].get_name();
        std::string as = std::to_string(a + 1);
        if (capMap.contains(an) && capMap.at(an).contains(an)) {
            double cself = capMap.at(an).at(an);
            if (std::isfinite(cself) && cself > 0) {
                capString += "Cself_" + as + " P" + as + "+ P" + as + "- " + to_string(cself, 18) + "\n";
                emittedAny = true;
            }
        }
        for (size_t b = a + 1; b < numWindings; ++b) {
            std::string bn = fd[b].get_name();
            if (capMap.contains(an) && capMap.at(an).contains(bn)) {
                double cinter = capMap.at(an).at(bn);
                if (std::isfinite(cinter) && cinter > 0) {
                    capString += "Cwind_" + as + "_" + std::to_string(b + 1) +
                                 " P" + as + "+ P" + std::to_string(b + 1) + "+ " +
                                 to_string(cinter, 18) + "\n";
                    emittedAny = true;
                }
            }
        }
    }
    return emittedAny ? capString : "";
}




// Strip leading/trailing whitespace (and CR) in place.
static void trim_inplace(std::string& s) {
    size_t start = 0;
    while (start < s.size() && (static_cast<unsigned char>(s[start]) <= 32)) ++start;
    size_t end = s.size();
    while (end > start && (static_cast<unsigned char>(s[end - 1]) <= 32)) --end;
    if (start > 0 || end < s.size()) {
        s = s.substr(start, end - start);
    }
}

// Decide whether a CSV line should be ignored entirely (blank or comment).
// Comment prefixes accepted: '#', '!', '*' (SPICE), "//".
static bool is_skippable_csv_line(const std::string& line) {
    if (line.empty()) return true;
    size_t i = 0;
    while (i < line.size() && static_cast<unsigned char>(line[i]) <= 32) ++i;
    if (i == line.size()) return true;
    char c = line[i];
    if (c == '#' || c == '!' || c == '*') return true;
    if (c == '/' && i + 1 < line.size() && line[i + 1] == '/') return true;
    return false;
}

// Strip a UTF-8 BOM (EF BB BF) if present at the start of the string.
static void strip_utf8_bom(std::string& s) {
    if (s.size() >= 3 &&
        static_cast<unsigned char>(s[0]) == 0xEF &&
        static_cast<unsigned char>(s[1]) == 0xBB &&
        static_cast<unsigned char>(s[2]) == 0xBF) {
        s.erase(0, 3);
    }
}

std::vector<std::string> CircuitSimulationReader::split_fields(const std::string& line, char separator) {
    // Mirrors getline(ss, token, separator) but treats a separator as literal
    // when it sits inside a double-quoted field or inside balanced parentheses.
    // Parenthesis-awareness is what keeps LTspice differential probes like
    // V(N009,d) from being split into two columns on a comma-separated file.
    std::vector<std::string> fields;
    std::string current;
    bool inQuotes = false;
    int parenDepth = 0;
    for (char c : line) {
        if (c == '"') {
            inQuotes = !inQuotes;
            current.push_back(c);  // keep quotes; the caller strips surrounding ones
            continue;
        }
        if (!inQuotes) {
            if (c == '(') {
                parenDepth++;
            }
            else if (c == ')' && parenDepth > 0) {
                parenDepth--;
            }
        }
        if (c == separator && !inQuotes && parenDepth == 0) {
            fields.push_back(current);
            current.clear();
            continue;
        }
        current.push_back(c);
    }
    fields.push_back(current);
    return fields;
}

void CircuitSimulationReader::process_line(std::string line, char separator) {
    process_line_with_context(line, separator, 0);
}

void CircuitSimulationReader::process_line_with_context(const std::string& line, char separator, size_t lineNumber) {
    if (_columns.size() == 0) {
        // Header row: collect column names. Strip CR and surrounding quotes from each token.
        for (std::string token : split_fields(line, separator)) {
            // Strip CR anywhere in the token (handles CRLF files).
            token.erase(std::remove(token.begin(), token.end(), '\r'), token.end());
            // Strip surrounding quotes only (preserves quoted commas/separators within names).
            if (token.size() >= 2 && token.front() == '"' && token.back() == '"') {
                token = token.substr(1, token.size() - 2);
            }
            // Trim whitespace; many exports pad headers with spaces.
            trim_inplace(token);
            if (token.empty()) continue;
            CircuitSimulationSignal circuitSimulationSignal;
            circuitSimulationSignal.name = token;
            _columns.push_back(circuitSimulationSignal);
        }
    }
    else {
        size_t currentColumnIndex = 0;
        for (std::string token : split_fields(line, separator)) {
            // Same CR/quote/whitespace cleanup as the header so a stray CR in the
            // last field of a CRLF file does not break std::stod.
            token.erase(std::remove(token.begin(), token.end(), '\r'), token.end());
            if (token.size() >= 2 && token.front() == '"' && token.back() == '"') {
                token = token.substr(1, token.size() - 2);
            }
            trim_inplace(token);
            if (token.empty()) {
                // An empty cell means "no sample for this column on this row".
                // The column index must still advance, otherwise every later
                // value on the row lands in the wrong column.
                currentColumnIndex++;
                continue;
            }
            if (currentColumnIndex >= _columns.size()) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "Data row " + std::to_string(lineNumber) +
                    " has more columns than the header (" +
                    std::to_string(_columns.size()) + ")");
            }
            double value;
            try {
                size_t consumed = 0;
                value = std::stod(token, &consumed);
                (void)consumed;
            }
            catch (const std::exception&) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "Could not parse number on line " + std::to_string(lineNumber) +
                    ", column " + std::to_string(currentColumnIndex + 1) + " (\"" +
                    _columns[currentColumnIndex].name + "\"): \"" + token + "\"");
            }
            _columns[currentColumnIndex].data.push_back(value);
            currentColumnIndex++;
        }
    }
}

CircuitSimulationReader::CircuitSimulationReader(std::string filePathOrFile, bool forceFile) {
    char separator = '\0';
    std::string line;
    size_t lineNumber = 0;

    std::filesystem::path path = filePathOrFile;

    auto handle_line = [&](std::string& rawLine) {
        ++lineNumber;
        // Strip BOM on the very first line; harmless if absent.
        if (lineNumber == 1) {
            strip_utf8_bom(rawLine);
        }
        if (is_skippable_csv_line(rawLine)) {
            return;
        }
        if (separator == '\0') {
            separator = guess_separator(rawLine);
        }
        process_line_with_context(rawLine, separator, lineNumber);
    };

    if (path.has_parent_path() && !forceFile) {
        std::ifstream is(filePathOrFile);
        if (!is.is_open()) {
            if (!std::filesystem::exists(filePathOrFile)) {
                throw InvalidInputException(ErrorCode::MISSING_DATA,
                    "File not found: " + filePathOrFile);
            }
            throw InvalidInputException(ErrorCode::MISSING_DATA,
                "Could not open file (check permissions): " + filePathOrFile);
        }
        while (getline(is, line)) {
            handle_line(line);
        }
        is.close();
    }
    else {
        std::stringstream ss(filePathOrFile);
        while (std::getline(ss, line, '\n')) {
            handle_line(line);
        }
    }

    if (_columns.empty()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "No columns found in input — the file may be empty, contain only comments, or use an unsupported separator");
    }

    _time = find_time(_columns);

    // Collapse repeated consecutive timestamps, keeping the LAST sample at each
    // timestamp. Simulator exports that round the time column to a few
    // significant figures emit several rows sharing one timestamp; without this
    // the time axis is not strictly increasing and period extraction would fail
    // ("no time column found"). Keeping the last value matches the most recent
    // state at that instant.
    {
        const std::vector<double>& times = _time.data;
        std::vector<size_t> keepIndexes;
        keepIndexes.reserve(times.size());
        for (size_t index = 0; index < times.size(); ++index) {
            if (index + 1 == times.size() || times[index] != times[index + 1]) {
                keepIndexes.push_back(index);
            }
        }
        if (keepIndexes.size() != times.size()) {
            for (auto& column : _columns) {
                if (column.data.size() != times.size()) {
                    continue;
                }
                std::vector<double> deduped;
                deduped.reserve(keepIndexes.size());
                for (size_t idx : keepIndexes) {
                    deduped.push_back(column.data[idx]);
                }
                column.data = std::move(deduped);
            }
            _time = find_time(_columns);
        }
    }
}

std::vector<std::string> CircuitSimulationReader::extract_column_names() {
    std::vector<std::string> names;
    for (auto column : _columns) {
        names.push_back(column.name);
    }

    return names;
}

bool CircuitSimulationReader::can_be_time(std::vector<double> data) {
    if (data.size() == 0) {
        throw std::invalid_argument("vector data cannot be empty");
    }
    if (data.size() == 1) {
        return false;
    }

    // Accept a monotonic non-decreasing axis. Time columns exported from
    // circuit simulators are frequently rounded to a few significant figures,
    // which yields repeated consecutive timestamps. Reject only a strictly
    // decreasing step (that is not a time axis), and require an overall
    // increase so a constant/zero column is never mistaken for time. Repeated
    // timestamps are collapsed after detection (keeping the last sample) so
    // downstream period extraction still sees a strictly increasing axis.
    for (size_t index = 1; index < data.size(); index++) {
        if (data[index - 1] > data[index]) {
            return false;
        }
    }

    return data.back() > data.front();
}

bool CircuitSimulationReader::can_be_voltage(std::vector<double> data, double limit) {
    if (data.size() == 0) {
        throw std::invalid_argument("vector data cannot be empty");
    }
    if (data.size() == 1) {
        return false;
    }

    std::vector<double> distinctValues;
    std::vector<size_t> distinctValuesCount;
    for (size_t index = 0; index < data.size(); index++) {
        bool isDistinct = true;
        for (size_t distincIndex = 0; distincIndex < distinctValues.size(); distincIndex++) {
            double distinctValue = distinctValues[distincIndex];
            double error = fabs(distinctValue - data[index]) / std::max(fabs(data[index]), fabs(distinctValue));
            if (std::isnan(error)) {
                error = fabs(distinctValue - data[index]);
            }
            if (error <= limit) {
                isDistinct = false;
                distinctValuesCount[distincIndex]++;
            }
        }
        if (isDistinct) {
            distinctValues.push_back(data[index]);
            distinctValuesCount.push_back(1);
        }
    }



    size_t significantDistinctValues = 0;
    for (size_t distincIndex = 0; distincIndex < distinctValuesCount.size(); distincIndex++) {
        if (distinctValuesCount[distincIndex] > data.size() * limit) {
            significantDistinctValues++;
        }
    }

    bool result = significantDistinctValues == 2 || significantDistinctValues == 3;
    return result;
}

bool CircuitSimulationReader::can_be_current(std::vector<double> data, double limit) {
    std::vector<double> diffValues;
    for (size_t index = 0; index < data.size(); index++) {
        double diff;
        if (index == 0) {
            diff = data[index] - data[data.size() - 1];
        }
        else {
            diff = data[index] - data[index - 1];
        }
        diffValues.push_back(diff);
    }

    return can_be_voltage(diffValues, limit);
}

char CircuitSimulationReader::guess_separator(std::string line){
    // Pick the candidate that produces the most columns (≥2), counting only
    // separator characters that appear *outside* double-quoted fields AND
    // outside balanced parentheses. This lets us correctly detect ';' as the
    // separator in
    //     "Voltage, primary";"Current, primary";"Time / s"
    // (commas live inside quotes) and ',' as the separator in an LTspice header
    //     time,V(N009,d),V(n016),I(L1)
    // where the comma inside the differential probe V(N009,d) must not be
    // counted — a naive std::count would over-count commas in both cases.
    static const std::vector<char> possibleSeparators = {',', ';', '\t', '|'};

    auto count_separators = [&](char sep) -> size_t {
        size_t count = 0;
        bool inQuotes = false;
        int parenDepth = 0;
        for (char c : line) {
            if (c == '"') {
                inQuotes = !inQuotes;
                continue;
            }
            if (!inQuotes) {
                if (c == '(') {
                    parenDepth++;
                }
                else if (c == ')' && parenDepth > 0) {
                    parenDepth--;
                }
            }
            if (!inQuotes && parenDepth == 0 && c == sep) {
                ++count;
            }
        }
        return count;
    };

    char bestSeparator = '\0';
    size_t bestColumns = 1;

    for (auto separator : possibleSeparators) {
        size_t numberColumns = count_separators(separator) + 1;
        if (numberColumns > bestColumns) {
            bestColumns = numberColumns;
            bestSeparator = separator;
        }
    }

    if (bestSeparator != '\0') {
        return bestSeparator;
    }

    throw InvalidInputException(ErrorCode::INVALID_INPUT,
        "No column separator found (expected one of comma, semicolon, tab, or pipe)");
}

Waveform CircuitSimulationReader::get_one_period(Waveform waveform, double frequency, bool sample, bool alignToZeroCrossing) {
    double period = 1.0 / frequency;
    if (!waveform.get_time()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Missing time data");
    }

    auto time = waveform.get_time().value();
    auto data = waveform.get_data();

    double periodEnd = time.back();
    double periodStart = periodEnd - period;
    size_t periodStartIndex = 0;
    size_t periodStopIndex = -1;

    if (_periodStartIndex && _periodStopIndex) {
        periodStartIndex = _periodStartIndex.value();
        periodStopIndex = _periodStopIndex.value();
    }
    else {
        for (int i = time.size() - 1; i >= 0; --i)
        {
            if (time[i] <= periodStart) {
                periodStartIndex = i;
                break;
            }
        }

        // Only search for zero crossing if alignToZeroCrossing is true
        if (alignToZeroCrossing) {
            double previousData = data[periodStartIndex];
            for (int i = periodStartIndex - 1; i >= 0; --i)
            {
                if ((data[i] >= 0 && previousData <= 0) || (data[i] <= 0 && previousData >= 0)) {
                    periodStartIndex = i;
                    periodStart = time[i];
                    break;
                }
                previousData = data[i];
            }
        }


        for (size_t i = periodStartIndex; i < time.size(); ++i)
        {
            if (time[i] >= periodStart + period) {
                periodStopIndex = i + 1;
                break;
            }
            else {

            }
        }
        // If the alignment back-walk pulled periodStart so far back that
        // periodStart + period > time.back(), the forward loop above never
        // assigns periodStopIndex (still SIZE_MAX from -1 init). Constructing
        // a vector with [begin+startIdx, begin+SIZE_MAX) trips length_error
        // deep in libstdc++. Surface a meaningful diagnostic at the source.
        if (periodStopIndex == static_cast<size_t>(-1)) {
            throw std::runtime_error(
                "CircuitSimulationReader::get_one_period: insufficient simulation "
                "data — could not find a complete period of " + std::to_string(period) +
                " s starting at t=" + std::to_string(periodStart) +
                " s (last sample at t=" + std::to_string(time.back()) +
                " s). Increase the converter's numSteadyStatePeriods or "
                "numPeriodsToExtract so the SPICE output covers at least one "
                "more switching period than the requested extraction window.");
        }
        _periodStartIndex = periodStartIndex;
        _periodStopIndex = periodStopIndex;
    }

    auto periodData = std::vector<double>(data.begin() + periodStartIndex, data.begin() + periodStopIndex);
    auto periodTime = std::vector<double>(time.begin() + periodStartIndex, time.begin() + periodStopIndex);


    double offset = periodTime[0];
    for (size_t index = 0; index < periodTime.size(); index++) {
        periodTime[index] -= offset;
    }



    Waveform newWaveform;
    newWaveform.set_data(periodData);
    newWaveform.set_time(periodTime);

    if (sample) {
        auto sampledWaveform = Inputs::calculate_sampled_waveform(newWaveform, frequency);
        return sampledWaveform;
    }
    else {
        return newWaveform;
    }
}

CircuitSimulationReader::CircuitSimulationSignal CircuitSimulationReader::find_time(std::vector<CircuitSimulationReader::CircuitSimulationSignal> columns) {
    for (auto column : columns) {
        if (can_be_time(column.data)) {
            column.type = DataType::TIME;
            return column;
        }
    }
    throw InvalidInputException(ErrorCode::MISSING_DATA, "no time column found");
}

Waveform CircuitSimulationReader::extract_waveform(CircuitSimulationReader::CircuitSimulationSignal signal, double frequency, bool sample) {
    Waveform waveform;
    waveform.set_data(signal.data);
    waveform.set_time(_time.data);
    return get_one_period(waveform, frequency, sample);
}

std::vector<int> get_numbers_in_string(std::string s) {
    std::regex regex(R"(\d+)");   // matches a sequence of digits
    std::vector<int> numbers;
    std::smatch match;
    while (std::regex_search(s, match, regex)) {
        numbers.push_back(std::stoi(match.str()));
        s = match.suffix();
    }

    return numbers;
}

bool CircuitSimulationReader::extract_winding_indexes(size_t numberWindings) {
    size_t numberFoundIndexes = 0;
    std::vector<size_t> indexes;
    // Insertion order matters: longer / more-specific labels must be checked
    // before short single-letter fallbacks, otherwise "I_trafo_HV" would match
    // "a" (from "trafo") before "HV". Using std::map iterated alphabetically,
    // which is what the previous code did by accident.
    const std::vector<std::pair<std::string, size_t>> windingLabels = {
        {"pri", 0},
        {"sec", 1},
        {"aux", 2},
        {"ter", 2},
        {"HV", 0},
        {"LV", 1},
        {"a", 0},
        {"b", 1},
        {"c", 2},
    };

    std::vector<CircuitSimulationSignal> columnsWithIndexes;
    std::vector<CircuitSimulationSignal> columnsWithResetIndexes;

    for (auto column : _columns) {
        column.windingIndex = -1;
        if (!can_be_time(column.data)) {
            auto numbersInColumnName = get_numbers_in_string(column.name);
            if (numbersInColumnName.size() > 0) {
                numberFoundIndexes++;
                indexes.push_back(numbersInColumnName.back());
                column.windingIndex = numbersInColumnName.back();
            }
            else {
                bool found = false;
                for (auto [label, windingIndex] : windingLabels) {
                    if (column.name.find(label) != std::string::npos)  {
                        numberFoundIndexes++;
                        indexes.push_back(windingIndex);
                        column.windingIndex = windingIndex;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    column.windingIndex = 0;
                    indexes.push_back(0);
                }
            }
        }
        columnsWithIndexes.push_back(column);
    }

    std::sort(indexes.begin(),indexes.end());
    indexes.resize(distance(indexes.begin(),std::unique(indexes.begin(),indexes.end())));

    for (size_t index = 0; index < indexes.size(); index++) {
        for (auto column : columnsWithIndexes) {
            if (column.windingIndex == indexes[index]) {
                column.windingIndex = index;
                columnsWithResetIndexes.push_back(column);
            }
        }
    }

    _columns = columnsWithResetIndexes;
    return numberFoundIndexes >= numberWindings;
}

std::vector<double> rolling_window_filter(std::vector<double> data) {

    size_t rollingFactorDividend = 192;
    size_t rollingFactor = data.size() / rollingFactorDividend;

    if (rollingFactor < 1){
        rollingFactor = 1;
    }

    for (size_t i = 0; i < rollingFactor - 1; ++i) {
        data.push_back(data[i]);
    }
    std::vector<double> ones(rollingFactor, 1.0);

    auto cleanData = convolution_valid(data, ones);
    for (size_t i = 0; i < cleanData.size(); ++i) {
        cleanData[i] /= rollingFactor;
    }
    return cleanData;
}

// Helper: check if a character is a word boundary (non-alphanumeric or start/end of string)
static bool is_word_boundary(const std::string& str, size_t pos) {
    if (pos == 0 || pos >= str.size()) return true;
    char c = str[pos - 1];
    return !std::isalnum(static_cast<unsigned char>(c)) && c != '_';
}

// Helper: check if alias matches at a word boundary in the name
static bool matches_at_word_boundary(const std::string& name, const std::string& prefix) {
    size_t pos = 0;
    while ((pos = name.find(prefix, pos)) != std::string::npos) {
        // Check that the match is at a word boundary (preceded by non-alnum or start of string)
        if (is_word_boundary(name, pos)) {
            // For single-char prefixes, also check that the next char (if any) is uppercase or non-alpha
            // This ensures "I" matches "Ipri" but not the "I" in "Inductor"
            size_t endPos = pos + prefix.size();
            if (endPos >= name.size()) return true;
            char nextChar = name[endPos];
            // If prefix is a single uppercase letter, the next char should not be lowercase
            // (to avoid matching "I" in "Inductor" or "V" in "Voltage")
            if (prefix.size() == 1 && std::isupper(static_cast<unsigned char>(prefix[0]))) {
                if (std::islower(static_cast<unsigned char>(nextChar))) {
                    pos++;
                    continue;
                }
            }
            return true;
        }
        pos++;
    }
    return false;
}

std::optional<CircuitSimulationReader::DataType> CircuitSimulationReader::guess_type_by_name(std::string name) {
    // 1. Check time aliases first (substring match is fine for multi-char aliases)
    for (const auto& timeAlias : _timeAliases) {
        if (name.find(timeAlias) != std::string::npos)  {
            return CircuitSimulationReader::DataType::TIME;
        }
    }

    // 2. Check magnetizing current BEFORE regular current (more specific first)
    for (const auto& magAlias : _magnetizingCurrentAliases) {
        if (name.find(magAlias) != std::string::npos)  {
            return CircuitSimulationReader::DataType::MAGNETIZING_CURRENT;
        }
    }

    // 3. Check multi-char current aliases (substring match)
    for (const auto& currentAlias : _currentAliases) {
        if (name.find(currentAlias) != std::string::npos)  {
            return CircuitSimulationReader::DataType::CURRENT;
        }
    }

    // 4. Check multi-char voltage aliases (substring match)
    for (const auto& voltageAlias : _voltageAliases) {
        if (name.find(voltageAlias) != std::string::npos)  {
            return CircuitSimulationReader::DataType::VOLTAGE;
        }
    }

    // 5. Check single-char prefixes with word-boundary awareness
    for (const auto& prefix : _currentPrefixes) {
        if (matches_at_word_boundary(name, prefix)) {
            return CircuitSimulationReader::DataType::CURRENT;
        }
    }

    for (const auto& prefix : _voltagePrefixes) {
        if (matches_at_word_boundary(name, prefix)) {
            return CircuitSimulationReader::DataType::VOLTAGE;
        }
    }

    return std::nullopt;
}

bool CircuitSimulationReader::extract_column_types(double frequency) {
    bool result = true;
    std::vector<std::map<std::string, std::string>> columnNameToSignalPerWinding;
    std::vector<CircuitSimulationSignal> columnsWithTypes;

    for (auto column : _columns) {
        if (can_be_time(column.data)) {
            column.type = DataType::TIME;
        }
        else {

            auto waveform = extract_waveform(column, frequency, false);
            auto data = waveform.get_data();
            int timeout = 100;

            auto guessedType = guess_type_by_name(column.name);

            if (guessedType) {
                column.type = guessedType.value();
            }
            else {
                while (timeout > 0) {
                    if (can_be_current(data)) {
                        column.type = DataType::CURRENT;
                        break;
                    }
                    else if (can_be_voltage(data)) {
                        column.type = DataType::VOLTAGE;
                        break;
                    }
                    else {
                        data = rolling_window_filter(data);
                    }
                    timeout--;
                } 
                if (timeout == 0) {
                    column.type = DataType::UNKNOWN;
                    // throw std::runtime_error("Unknown type");
                }
            }
            columnsWithTypes.push_back(column);
        }
    }
    _columns = columnsWithTypes;
    return result;
}

std::vector<std::map<std::string, std::string>> CircuitSimulationReader::extract_map_column_names(size_t numberWindings, double frequency) {
    OperatingPoint operatingPoint;
    std::vector<std::map<std::string, std::string>> columnNameToSignalPerWinding;

    std::vector<OperatingPointExcitation> excitationsPerWinding;

    extract_winding_indexes(numberWindings);
    extract_column_types(frequency);
    for (size_t windingIndex = 0; windingIndex < numberWindings; windingIndex++) {
        std::map<std::string, std::string> columnNameToSignal;
        columnNameToSignal["time"] = _time.name;
        columnNameToSignal["current"] = "";
        columnNameToSignal["voltage"] = "";
        for (auto column : _columns) {
            if (column.windingIndex == windingIndex && column.type == DataType::CURRENT) {
                columnNameToSignal["current"] = column.name;
            }
            if (column.windingIndex == windingIndex && column.type == DataType::MAGNETIZING_CURRENT) {
                columnNameToSignal["magnetizingCurrent"] = column.name;
            }
            if (column.windingIndex == windingIndex && column.type == DataType::VOLTAGE) {
                columnNameToSignal["voltage"] = column.name;
            }
        }

        columnNameToSignalPerWinding.push_back(columnNameToSignal);
    }

    return columnNameToSignalPerWinding;
}

bool CircuitSimulationReader::assign_column_names(std::vector<std::map<std::string, std::string>> columnNames) {
    std::vector<CircuitSimulationSignal> assignedColumns;
    for (size_t windingIndex = 0; windingIndex < columnNames.size(); windingIndex++) {
        for (auto [columnType, columnName] : columnNames[windingIndex]) {

            for (auto column : _columns) {
                if (column.name == columnName && columnType == "current") {
                    column.type = DataType::CURRENT;
                    column.windingIndex = windingIndex;
                    assignedColumns.push_back(column);
                }
                if (column.name == columnName && columnType == "magnetizingCurrent") {
                    column.type = DataType::MAGNETIZING_CURRENT;
                    column.windingIndex = windingIndex;
                    assignedColumns.push_back(column);
                }
                if (column.name == columnName && columnType == "voltage") {
                    column.type = DataType::VOLTAGE;
                    column.windingIndex = windingIndex;
                    assignedColumns.push_back(column);
                }
                if (column.name == columnName && columnType == "time") {
                    column.type = DataType::TIME;
                    column.windingIndex = windingIndex;
                    assignedColumns.push_back(column);
                }
            }
        }
    }

    _columns = assignedColumns;

    return true;
}

OperatingPoint CircuitSimulationReader::extract_operating_point(size_t numberWindings, double frequency, std::optional<std::vector<std::map<std::string, std::string>>> mapColumnNames, double ambientTemperature) {
    OperatingPoint operatingPoint;

    std::vector<OperatingPointExcitation> excitationsPerWinding;
    if (!mapColumnNames) {
        extract_winding_indexes(numberWindings);
        extract_column_types(frequency);
    }
    else {
        assign_column_names(mapColumnNames.value());
    }

    for (size_t windingIndex = 0; windingIndex < numberWindings; windingIndex++) {
        OperatingPointExcitation excitation;
        for (auto column : _columns) {
            excitation.set_frequency(frequency);
            if (column.windingIndex == windingIndex && column.type == DataType::CURRENT) {
                auto waveform = extract_waveform(column, frequency);
                SignalDescriptor current;
                current.set_waveform(waveform);
                // Calculate processed data for current
                auto currentProcessed = Inputs::calculate_processed_data(waveform, frequency, true);
                current.set_processed(currentProcessed);
                // Harmonics: downstream consumers (Inputs::prune_harmonics, core/copper
                // loss models) require .harmonics() to be set on every SignalDescriptor
                // returned from a simulation. Skipping this caused bad_optional_access
                // in DAB / resonant simulated operating-point flows.
                current.set_harmonics(Inputs::calculate_harmonics_data(waveform, frequency));
                excitation.set_current(current);
            }
            if (column.windingIndex == windingIndex && column.type == DataType::MAGNETIZING_CURRENT) {
                auto waveform = extract_waveform(column, frequency);
                SignalDescriptor current;
                current.set_waveform(waveform);
                // Calculate processed data for magnetizing current
                auto currentProcessed = Inputs::calculate_processed_data(waveform, frequency, true);
                current.set_processed(currentProcessed);
                current.set_harmonics(Inputs::calculate_harmonics_data(waveform, frequency));
                excitation.set_magnetizing_current(current);
            }
            if (column.windingIndex == windingIndex && column.type == DataType::VOLTAGE) {
                auto waveform = extract_waveform(column, frequency);
                SignalDescriptor voltage;
                voltage.set_waveform(waveform);
                // Calculate processed data for voltage
                auto voltageProcessed = Inputs::calculate_processed_data(waveform, frequency, true);
                voltage.set_processed(voltageProcessed);
                voltage.set_harmonics(Inputs::calculate_harmonics_data(waveform, frequency));
                excitation.set_voltage(voltage);
            }
        }
        excitationsPerWinding.push_back(excitation);
    }

    // Remove trailing empty excitations (no current, no voltage, no magnetizing current)
    // This handles the case where the user requests more windings than the data supports
    while (excitationsPerWinding.size() > 1) {
        auto& last = excitationsPerWinding.back();
        if (!last.get_current() && !last.get_voltage() && !last.get_magnetizing_current()) {
            excitationsPerWinding.pop_back();
        }
        else {
            break;
        }
    }

    operatingPoint.set_excitations_per_winding(excitationsPerWinding);
    OperatingConditions conditions;
    conditions.set_cooling(std::nullopt);    
    conditions.set_ambient_temperature(ambientTemperature);
    operatingPoint.set_conditions(conditions);

    return operatingPoint;
}




// =============================================================================
// FRACPOLE SKIN-EFFECT MODEL
// =============================================================================

double CircuitSimulatorExporter::fracpole_skin_model(double x[], double frequency) {
    double h = x[0];
    double alpha = x[1];
    if (h < 0) return 1e30;
    if (alpha < 0.01) return 1e30;
    if (alpha > 1.5) return 1e30;
    return h * pow(2 * std::numbers::pi * frequency, alpha);
}

void CircuitSimulatorExporter::fracpole_skin_func(double *p, double *x, int m, int n, void *data) {
    double* dcResistanceAndFrequencies = static_cast<double*>(data);
    double dcResistance = dcResistanceAndFrequencies[0];
    for (int i = 0; i < n; ++i) {
        x[i] = dcResistance + fracpole_skin_model(p, dcResistanceAndFrequencies[i + 1]);
    }
}

double CircuitSimulatorExporter::fracpole_model(double x[], int numCoeffs, double frequency, double dcResistance) {
    double w = 2 * std::numbers::pi * frequency;
    std::complex<double> Z(dcResistance, 0);
    int maxStages = numCoeffs / 2;
    for (int k = 0; k < maxStages; k++) {
        if (x[2*k] <= 0 || x[2*k+1] <= 0) break;
        auto R_k = std::complex<double>(x[2*k], 0);
        auto jwL_k = std::complex<double>(0, w * x[2*k+1]);
        Z += parallel(jwL_k, R_k);
    }
    return Z.real();
}

CircuitSimulatorExporterCurveFittingModes CircuitSimulatorExporter::resolve_curve_fitting_mode(
    CircuitSimulatorExporterCurveFittingModes mode) {
    if (mode == CircuitSimulatorExporterCurveFittingModes::AUTO) {
        auto& settings = OpenMagnetics::Settings::GetInstance();
        return static_cast<CircuitSimulatorExporterCurveFittingModes>(settings.get_circuit_simulator_curve_fitting_mode());
    }
    return mode;
}

std::vector<FractionalPoleNetwork> CircuitSimulatorExporter::calculate_fracpole_networks_per_winding(
    Magnetic magnetic, double temperature, FractionalPoleOptions opts) {

    auto& settings = OpenMagnetics::Settings::GetInstance();
    auto settingsOpts = settings.get_circuit_simulator_fracpole_options();
    if (opts.f0 <= 0) opts.f0 = settingsOpts ? (*settingsOpts)["f0"] : 0.1;
    if (opts.f1 <= opts.f0) opts.f1 = settingsOpts ? (*settingsOpts)["f1"] : 1e7;
    if (opts.lumpsPerDecade <= 0) opts.lumpsPerDecade = settingsOpts ? (*settingsOpts)["lumpsPerDecade"] : 1.5;
    opts.alpha = 0.5;  // Fixed by skin effect physics
    opts.profile = FracpoleProfile::DD;  // Recommended for alpha <= 0.5

    auto coil = magnetic.get_coil();
    std::vector<FractionalPoleNetwork> networksPerWinding;

    for (size_t windingIndex = 0;
         windingIndex < coil.get_functional_description().size();
         ++windingIndex) {
        double dcResistance = WindingLosses::calculate_effective_resistance_of_winding(
            magnetic, windingIndex, opts.f0, temperature);

        opts.coef = FractionalPole::fit_coef_from_reference(
            opts.f0, dcResistance, opts.alpha, opts.f0, opts.f1,
            opts.lumpsPerDecade, opts.profile);

        FractionalPoleNetwork net = FractionalPole::generate(opts);
        networksPerWinding.push_back(net);
    }
    return networksPerWinding;
}

FractionalPoleNetwork CircuitSimulatorExporter::calculate_core_fracpole_network(
    Magnetic magnetic, double temperature, FractionalPoleOptions opts) {

    auto& settings = OpenMagnetics::Settings::GetInstance();
    auto settingsOpts = settings.get_circuit_simulator_fracpole_options();
    if (opts.f0 <= 0) opts.f0 = settingsOpts ? (*settingsOpts)["f0"] : 0.1;
    if (opts.f1 <= opts.f0) opts.f1 = settingsOpts ? (*settingsOpts)["f1"] : 1e7;
    if (opts.lumpsPerDecade <= 0) opts.lumpsPerDecade = settingsOpts ? (*settingsOpts)["lumpsPerDecade"] : 1.5;

    // Extract Steinmetz alpha from core material. A material without Steinmetz
    // data cannot be silently modeled as "typical ferrite" — let it throw.
    auto core = magnetic.get_core();
    auto material = core.resolve_material();
    double refFreq = std::sqrt(opts.f0 * opts.f1);
    auto steinmetzDatum = CoreLossesModel::get_steinmetz_coefficients(material, refFreq);
    double steinmetzAlpha = steinmetzDatum.get_alpha();

    // alpha for fracpole: core resistance R_core ∝ f^(steinmetz_alpha - 1)
    opts.alpha = std::clamp(steinmetzAlpha - 1.0, 0.05, 0.95);

    // Choose profile: DD for alpha <= 0.5, FF for alpha > 0.5 (per Kundert recommendation)
    opts.profile = (opts.alpha <= 0.5) ? FracpoleProfile::DD : FracpoleProfile::FF;

    Curve2D coreResData = Sweeper().sweep_core_resistance_over_frequency(
        magnetic, opts.f0, opts.f1, 5, temperature);
    auto coreResVec = coreResData.get_y_points();
    auto freqVec    = coreResData.get_x_points();

    if (coreResVec.empty()) {
        throw std::runtime_error("calculate_core_fracpole_network: core resistance sweep returned no data, cannot fit the fracpole network");
    }

    size_t midIdx = coreResVec.size() / 2;
    double R_ref = coreResVec[midIdx];
    double f_ref_actual = freqVec[midIdx];

    opts.coef = FractionalPole::fit_coef_from_reference(
        f_ref_actual, R_ref, opts.alpha, opts.f0, opts.f1,
        opts.lumpsPerDecade, opts.profile);

    return FractionalPole::generate(opts);
}


GseCoreLossParams CircuitSimulatorExporter::calculate_gse_core_loss_params(
    Magnetic magnetic, double frequency, double temperature) {

    GseCoreLossParams p;

    auto core = magnetic.get_core();
    auto material = core.resolve_material();
    // Steinmetz k / alpha / beta for the material at this frequency band. A material without
    // Steinmetz data cannot supply a large-signal core-loss model -> report invalid and let the
    // caller keep the linear small-signal ladder.
    SteinmetzCoreLossesMethodRangeDatum steinmetzDatum;
    try {
        steinmetzDatum = CoreLossesModel::get_steinmetz_coefficients(material, frequency);
    } catch (...) {
        return p;
    }
    double alpha = steinmetzDatum.get_alpha();
    double beta = steinmetzDatum.get_beta();
    double k = steinmetzDatum.get_k();
    if (!(std::isfinite(alpha) && std::isfinite(beta) && std::isfinite(k) &&
          alpha > 0 && beta > 0 && k > 0 && beta >= alpha)) {
        return p;  // beta < alpha would make |B|^(beta-alpha) singular at B=0; refuse it
    }

    // GSE coefficient k1: calibrated so that for a sinusoidal B the cycle-average of
    // p(t) = k1*|dB/dt|^alpha*|B|^(beta-alpha) equals the Steinmetz loss k*f^alpha*B_peak^beta.
    //   k1 = k / [ (2*pi)^(alpha-1) * integral_0^{2*pi} |cos t|^alpha * |sin t|^(beta-alpha) dt ]
    // (identical construction to CoreLossesIGSEModel::get_ki, but with |sin t|^(beta-alpha) in the
    // integrand instead of the constant 2^(beta-alpha): GSE uses the instantaneous |B|, iGSE the
    // per-cycle peak-to-peak swing.)
    auto& settings = Settings::GetInstance();
    size_t nPoints = settings.get_inputs_number_points_sampled_waveforms();
    if (nPoints < 8) nPoints = 128;
    double integral = 0.0;
    for (size_t i = 0; i < nPoints; ++i) {
        double theta = static_cast<double>(i) * 2.0 * std::numbers::pi / static_cast<double>(nPoints);
        integral += std::pow(std::fabs(std::cos(theta)), alpha) *
                    std::pow(std::fabs(std::sin(theta)), beta - alpha) *
                    2.0 * std::numbers::pi / static_cast<double>(nPoints);
    }
    if (!(integral > 0.0)) {
        return p;
    }
    double k1 = k / (std::pow(2.0 * std::numbers::pi, alpha - 1.0) * integral);

    // Core geometry: effective volume + effective area, and the main-winding turns.
    auto processedDescription = core.get_processed_description();
    if (!processedDescription) {
        return p;
    }
    auto effectiveParameters = processedDescription.value().get_effective_parameters();
    double A_e = effectiveParameters.get_effective_area();
    double V_e = effectiveParameters.get_effective_volume();
    auto coil = magnetic.get_coil();
    if (coil.get_functional_description().empty()) {
        return p;
    }
    double N = static_cast<double>(coil.get_functional_description()[0].get_number_turns());
    if (!(A_e > 0 && V_e > 0 && N > 0)) {
        return p;
    }

    // I_loss = P/V_L with P = V_e*k1*|dB/dt|^alpha*|B|^(beta-alpha), dB/dt = V_L/(N*A_e),
    // B = lambda/(N*A_e):  I_loss = [V_e*k1/(N*A_e)^beta] * sgn(V_L)*|V_L|^(alpha-1)*|lambda|^(beta-alpha).
    p.gc = V_e * k1 / std::pow(N * A_e, beta);
    p.alpha = alpha;
    p.beta = beta;
    p.nTurns = N;
    p.effectiveArea = A_e;
    p.valid = std::isfinite(p.gc) && p.gc > 0;
    return p;
}


std::string emit_gse_core_loss_spice(
    const GseCoreLossParams& params,
    const std::string& windingIndex,
    const std::string& nodeIn,
    const std::string& nodeOut,
    bool quoteExpression) {

    if (!params.valid) {
        return "";
    }
    std::string fluxNode = "Node_Bflux_" + windingIndex;
    // Smooth (Lipschitz) realization of I = gc*sgn(V_L)*|V_L|^(alpha-1)*|lambda|^(beta-alpha),
    // written as V_L*(V_L^2+epsV^2)^((alpha-2)/2) * (lambda^2+epsL^2)^((beta-alpha)/2). The raw
    // |V_L|^(alpha-1) has an infinite slope at V_L=0 (every voltage zero-crossing) and stalls the
    // transient ("timestep too small"); the +eps^2 forms are bounded everywhere with negligible
    // bias (eps << any real winding voltage / flux linkage), so the element converges at full speed.
    double halfAlphaM2 = (params.alpha - 2.0) / 2.0;
    double halfBetaMAlpha = (params.beta - params.alpha) / 2.0;
    const double epsV2 = 1e-6;   // (1 mV)^2 voltage floor
    const double epsL2 = 1e-18;  // (1 nWb)^2 flux-linkage floor
    std::string vL = "V(" + nodeIn + "," + nodeOut + ")";
    std::string expr = to_string(params.gc, 15) + "*" + vL +
        "*pow(" + vL + "*" + vL + "+" + to_string(epsV2, 12) + "," + to_string(halfAlphaM2, 9) + ")" +
        "*pow(V(" + fluxNode + ")*V(" + fluxNode + ")+" + to_string(epsL2, 21) + "," +
        to_string(halfBetaMAlpha, 9) + ")";
    std::string q = quoteExpression ? "'" : "";

    std::string s;
    s += "* Behavioural GSE core loss (winding " + windingIndex +
         "): p=k1*|dB/dt|^alpha*|B|^(beta-alpha)\n";
    s += "* alpha=" + to_string(params.alpha, 4) + " beta=" + to_string(params.beta, 4) +
         " gc=" + to_string(params.gc, 15) + "\n";
    // Flux linkage lambda = integral(V_L) dt : winding voltage integrated into a 1 F cap node.
    s += "Gcl_flux_" + windingIndex + " 0 " + fluxNode + " " + nodeIn + " " + nodeOut + " 1\n";
    s += "Ccl_flux_" + windingIndex + " " + fluxNode + " 0 1\n";
    s += "Rcl_flux_" + windingIndex + " " + fluxNode + " 0 1e12\n";  // keep the flux node non-floating
    s += "Bcloss_" + windingIndex + " " + nodeIn + " " + nodeOut + " I=" + q + expr + q + "\n";
    return s;
}


} // namespace OpenMagnetics

