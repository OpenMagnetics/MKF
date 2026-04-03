#include "support/Utils.h"
#include "Defaults.h"
#include "processors/CircuitSimulatorInterface.h"
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

namespace OpenMagnetics {

// Saturation model parameters for circuit simulation
struct SaturationParameters {
    double Bsat;           // Saturation flux density [T]
    double mu_r;           // Relative permeability (initial)
    double Ae;             // Effective core area [m²]
    double le;             // Effective magnetic path length [m]
    double primaryTurns;   // Number of primary turns
    double Lmag;           // Magnetizing inductance [H]
    bool valid;            // Whether parameters are valid
    
    // Derived parameters for circuit models
    double Isat() const {
        // Current at which core saturates: I = H_sat * le / N
        // H_sat = B_sat / (mu_0 * mu_r)
        const double mu0 = 4e-7 * M_PI;
        double Hsat = Bsat / (mu0 * mu_r);
        return Hsat * le / primaryTurns;
    }
    
    double fluxLinkageSat() const {
        // Lambda_sat = N * Ae * Bsat
        return primaryTurns * Ae * Bsat;
    }
};

// Extract saturation parameters from a magnetic component
static SaturationParameters get_saturation_parameters(const Magnetic& magnetic, double temperature) {
    SaturationParameters params;
    params.valid = false;
    
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
        
        // Validate
        if (params.Bsat > 0 && params.mu_r > 0 && params.Ae > 0 && 
            params.le > 0 && params.primaryTurns > 0 && params.Lmag > 0) {
            params.valid = true;
        }
    } catch (...) {
        params.valid = false;
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

std::vector<double> CircuitSimulatorExporter::calculate_core_resistance_coefficients(Magnetic magnetic, double temperature) {
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
        dcResistanceAndfrequencies[0] = coreResistances[0];
        // coefficients[1] = 1e-7;
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

CircuitSimulatorExporterSimbaModel::CircuitSimulatorExporterSimbaModel() {
    // FIXED: Was seeding with random_device then immediately overwriting with time(0)
    // Now use random_device only, which provides proper entropy
    std::random_device rd;
    _gen = std::mt19937(rd());
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
    else
        throw ModelNotAvailableException("Unknown Circuit Simulator program, available options are: {SIMBA, NGSPICE, LTSPICE, PLECS}");
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

std::string CircuitSimulatorExporterSimbaModel::generate_id() {
    // generator for hex numbers from 0 to F
    std::uniform_int_distribution<int> dis(0, 15);
    
    const char hexChars[] = { '0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F' };
    const std::size_t hexCharsSize = sizeof(hexChars) / sizeof(hexChars[0]);
    std::string id = "";
    for (std::size_t i = 0; i < 8; ++i) {
        id += hexChars[static_cast<std::size_t>(dis(_gen)) % hexCharsSize];
    }
    id += "-";
    for (std::size_t i = 0; i < 4; ++i) {
        id += hexChars[static_cast<std::size_t>(dis(_gen)) % hexCharsSize];
    }
    id += "-";
    for (std::size_t i = 0; i < 4; ++i) {
        id += hexChars[static_cast<std::size_t>(dis(_gen)) % hexCharsSize];
    }
    id += "-";
    for (std::size_t i = 0; i < 4; ++i) {
        id += hexChars[static_cast<std::size_t>(dis(_gen)) % hexCharsSize];
    }
    id += "-";
    for (std::size_t i = 0; i < 12; ++i) {
        id += hexChars[static_cast<std::size_t>(dis(_gen)) % hexCharsSize];
    }
    return id;
}

ordered_json CircuitSimulatorExporterSimbaModel::create_device(std::string libraryName, std::vector<int> coordinates, int angle, std::string name) {
    ordered_json deviceJson;

    deviceJson["LibraryName"] = libraryName;
    deviceJson["Top"] = coordinates[1];
    deviceJson["Left"] = coordinates[0];
    deviceJson["Angle"] = angle;
    deviceJson["HF"] = false;
    deviceJson["VF"] = false;
    deviceJson["Disabled"] = false;
    deviceJson["Name"] = name;
    deviceJson["ID"] = generate_id();
    deviceJson["Parameters"]["Name"] = name;
    deviceJson["EnabledScopes"] = ordered_json::array();

    return deviceJson;
}

ordered_json CircuitSimulatorExporterSimbaModel::create_air_gap(std::vector<int> coordinates, double area, double length, int angle, std::string name) {
    ordered_json deviceJson = create_device("Air Gap", coordinates, angle, name);

    deviceJson["Parameters"]["RelativePermeability"] = "1";
    deviceJson["Parameters"]["Area"] = to_string(area, _precision);
    deviceJson["Parameters"]["Length"] = to_string(length, _precision);
    return deviceJson;
}
ordered_json CircuitSimulatorExporterSimbaModel::create_core(double initialPermeability, std::vector<int> coordinates, double area, double length, int angle, std::string name) {
    ordered_json deviceJson = create_device("Linear Core", coordinates, angle, name);

    deviceJson["Parameters"]["RelativePermeability"] = to_string(initialPermeability, _precision);
    deviceJson["Parameters"]["Area"] = to_string(area, _precision);
    deviceJson["Parameters"]["Length"] = to_string(length, _precision);
    return deviceJson;
}

ordered_json CircuitSimulatorExporterSimbaModel::create_winding(size_t numberTurns, std::vector<int> coordinates, int angle, std::string name) {
    ordered_json deviceJson = create_device("Winding", coordinates, angle, name);

    deviceJson["Parameters"]["NumberOfTurns"] = to_string(numberTurns, _precision);
    return deviceJson;
}

ordered_json CircuitSimulatorExporterSimbaModel::create_resistor(double resistance, std::vector<int> coordinates, int angle, std::string name) {

    ordered_json deviceJson = create_device("Resistor", coordinates, angle, name);

    deviceJson["Parameters"]["Value"] = to_string(resistance, _precision);
    return deviceJson;
}

ordered_json CircuitSimulatorExporterSimbaModel::create_inductor(double inductance, std::vector<int> coordinates, int angle, std::string name) {

    ordered_json deviceJson = create_device("Inductor", coordinates, angle, name);

    deviceJson["Parameters"]["Value"] = to_string(inductance, _precision);
    deviceJson["Parameters"]["Iinit"] = "0";
    return deviceJson;
}

ordered_json CircuitSimulatorExporterSimbaModel::create_capacitor(double capacitance, std::vector<int> coordinates, int angle, std::string name) {
    ordered_json deviceJson = create_device("Capacitor", coordinates, angle, name);
    deviceJson["Parameters"]["Value"] = to_string(capacitance, _precision);
    deviceJson["Parameters"]["Vinit"] = "0";
    return deviceJson;
}

std::pair<std::vector<ordered_json>, std::vector<ordered_json>> CircuitSimulatorExporterSimbaModel::create_ladder(std::vector<double> ladderCoefficients, std::vector<int> coordinates, std::string name) {
    std::vector<ordered_json> ladderJsons;
    std::vector<ordered_json> ladderConnectorsJsons;
    coordinates[0] -= 6;
    for (size_t ladderIndex = 0; ladderIndex < ladderCoefficients.size(); ladderIndex+=2) {
        std::string ladderIndexs = std::to_string(ladderIndex);
        auto ladderInductorJson = create_inductor(ladderCoefficients[ladderIndex + 1], coordinates, 180, name + " Ladder element " + std::to_string(ladderIndex));
        ladderJsons.push_back(ladderInductorJson);
        coordinates[1] -= 3;
        coordinates[0] += 3;
        auto ladderResistorJson = create_resistor(ladderCoefficients[ladderIndex], coordinates, 90, name + " Ladder element " + std::to_string(ladderIndex));
        ladderJsons.push_back(ladderResistorJson);
        coordinates[1] -= 3;
        coordinates[0] -= 3;
        {
            std::vector<int> finalConnectorTopCoordinates = {coordinates[0], coordinates[1] + 1};
            std::vector<int> finalConnectorBottomCoordinates = {coordinates[0], coordinates[1] + 7};
            auto ladderConnectorJson = create_connector(finalConnectorTopCoordinates, finalConnectorBottomCoordinates, "Connector between winding " + name + " ladder elements" + std::to_string(ladderIndex) + " and " + std::to_string(ladderIndex + 2));
            ladderConnectorsJsons.push_back(ladderConnectorJson);
        }
        if (ladderIndex == ladderCoefficients.size() - 2) {
            std::vector<int> finalConnectorTopCoordinates = {coordinates[0], coordinates[1] + 1};
            std::vector<int> finalConnectorBottomCoordinates = {coordinates[0] + 6, coordinates[1] + 1};
            auto ladderConnectorJson = create_connector(finalConnectorTopCoordinates, finalConnectorBottomCoordinates, "Connector between winding " + name + " ladder elements" + std::to_string(ladderIndex) + " and " + std::to_string(ladderIndex + 2));
            ladderConnectorsJsons.push_back(ladderConnectorJson);
        }
    }

    return {ladderJsons, ladderConnectorsJsons};
}

std::pair<std::vector<ordered_json>, std::vector<ordered_json>> CircuitSimulatorExporterSimbaModel::create_fracpole_ladder(
    const FractionalPoleNetwork& net, std::vector<int> coordinates, const std::string& name) {
    std::vector<ordered_json> ladderJsons;
    std::vector<ordered_json> ladderConnectorsJsons;
    coordinates[0] -= 6;

    for (size_t k = 0; k < net.stages.size(); ++k) {
        std::string ks = std::to_string(k);
        // Resistor (series, horizontal)
        auto resistorJson = create_resistor(net.stages[k].R * net.opts.coef, coordinates, 180,
                                            name + " Fracpole R" + ks);
        ladderJsons.push_back(resistorJson);
        coordinates[1] -= 3;
        coordinates[0] += 3;
        // Capacitor (shunt, vertical) - replaces the inductor in ladder mode
        auto capacitorJson = create_capacitor(net.stages[k].C / net.opts.coef, coordinates, 90,
                                              name + " Fracpole C" + ks);
        ladderJsons.push_back(capacitorJson);
        coordinates[1] -= 3;
        coordinates[0] -= 3;
        // Connector between stages
        {
            std::vector<int> top = {coordinates[0], coordinates[1] + 1};
            std::vector<int> bot = {coordinates[0], coordinates[1] + 7};
            auto connJson = create_connector(top, bot,
                "Fracpole connector " + name + " stage " + ks);
            ladderConnectorsJsons.push_back(connJson);
        }
        if (k == net.stages.size() - 1) {
            std::vector<int> top = {coordinates[0], coordinates[1] + 1};
            std::vector<int> bot = {coordinates[0] + 6, coordinates[1] + 1};
            auto connJson = create_connector(top, bot, "Fracpole final connector " + name);
            ladderConnectorsJsons.push_back(connJson);
        }
    }
    return {ladderJsons, ladderConnectorsJsons};
}

ordered_json CircuitSimulatorExporterSimbaModel::create_pin(std::vector<int> coordinates, int angle, std::string name) {
    return create_device("Electrical Pin", coordinates, angle, name);
}

ordered_json CircuitSimulatorExporterSimbaModel::create_magnetic_ground(std::vector<int> coordinates, int angle, std::string name) {
    return create_device("Magnetic Ground", coordinates, angle, name);
}

ordered_json CircuitSimulatorExporterSimbaModel::create_connector(std::vector<int> startingCoordinates, std::vector<int> endingCoordinates, std::string name) {
    ordered_json connectorJson;

    connectorJson["Name"] = name;
    connectorJson["Segments"] = ordered_json::array();

    if (startingCoordinates[0] == endingCoordinates[0] || startingCoordinates[1] == endingCoordinates[1]) {
        ordered_json segment;
        segment["StartX"] = startingCoordinates[0];
        segment["StartY"] = startingCoordinates[1];
        segment["EndX"] = endingCoordinates[0];
        segment["EndY"] = endingCoordinates[1];
        connectorJson["Segments"].push_back(segment);
    }
    else {
        {
            ordered_json segment;
            segment["StartX"] = startingCoordinates[0];
            segment["StartY"] = startingCoordinates[1];
            segment["EndX"] = endingCoordinates[0];
            segment["EndY"] = startingCoordinates[1];
            connectorJson["Segments"].push_back(segment);
        }
        {
            ordered_json segment;
            segment["StartX"] = endingCoordinates[0];
            segment["StartY"] = startingCoordinates[1];
            segment["EndX"] = endingCoordinates[0];
            segment["EndY"] = endingCoordinates[1];
            connectorJson["Segments"].push_back(segment);
        }
    }

    return connectorJson;
}

ordered_json CircuitSimulatorExporterSimbaModel::merge_connectors(ordered_json connectors) {
    auto mergeConnectors = connectors;
    bool merged = true;
    while (merged) {
        merged = false;
        for (size_t firstIndex = 0; firstIndex < mergeConnectors.size(); ++firstIndex) {
            bool sharedPoint = false;
            ordered_json mergedConnector;
            for (size_t secondIndex = firstIndex + 1; secondIndex < mergeConnectors.size(); ++secondIndex) {
                for (auto firstConnectorSegment : mergeConnectors[firstIndex]["Segments"]) {
                    for (auto secondConnectorSegment : mergeConnectors[secondIndex]["Segments"]) {
                        if ((firstConnectorSegment["StartX"] == secondConnectorSegment["StartX"] &&
                            firstConnectorSegment["StartY"] == secondConnectorSegment["StartY"]) ||
                            (firstConnectorSegment["EndX"] == secondConnectorSegment["EndX"] &&
                            firstConnectorSegment["EndY"] == secondConnectorSegment["EndY"])) {

                            sharedPoint = true;
                            break;
                        }
                    }
                    if (sharedPoint) {
                        break;
                    }
                }

                if (sharedPoint) {
                    mergedConnector["Segments"] = ordered_json::array();
                    mergedConnector["Name"] = "Merge of connector: " + std::string{mergeConnectors[firstIndex]["Name"]} + " with " + std::string{mergeConnectors[secondIndex]["Name"]};
                    for (auto firstConnectorSegment : mergeConnectors[firstIndex]["Segments"]) {
                        mergedConnector["Segments"].push_back(firstConnectorSegment);
                    }
                    for (auto secondConnectorSegment : mergeConnectors[secondIndex]["Segments"]) {
                        mergedConnector["Segments"].push_back(secondConnectorSegment);
                    }
                    mergeConnectors.erase(secondIndex);
                    mergeConnectors.erase(firstIndex);
                    mergeConnectors.push_back(mergedConnector);
                    merged = true;
                    break;
                }
            }

            if (sharedPoint) {
                break;
            }
        }
    }
    return mergeConnectors;
}



std::string CircuitSimulatorExporterSimbaModel::export_magnetic_as_subcircuit(Magnetic magnetic, double frequency, double temperature, std::optional<std::string> filePathOrFile, CircuitSimulatorExporterCurveFittingModes mode) {
    ordered_json simulation;
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();
    _scale = _modelSize / core.get_width();


    if (filePathOrFile) {
        try {
            if(!std::filesystem::exists(filePathOrFile.value())) {
                throw InvalidInputException(ErrorCode::MISSING_DATA, "File not found");
            }
            std::ifstream json_file(filePathOrFile.value());
            if(json_file.is_open()) {
                simulation = ordered_json::parse(json_file);
            }
            else {
                throw InvalidInputException(ErrorCode::MISSING_DATA, "File not found");
            }
        }
        catch(const InvalidInputException& re) {
            std::stringstream json_file(filePathOrFile.value());
            simulation = ordered_json::parse(json_file);
        }
    }

    if (!simulation.contains("Libraries")) {
        simulation["Libraries"] = ordered_json::array();
    }
    if (!simulation.contains("Designs")) {
        simulation["Designs"] = ordered_json::array();
    }
    ordered_json library;
    ordered_json device;
    library["LibraryItemName"] = "OpenMagnetics components";
    // library["Devices"] = ordered_json::array();

    device["LibraryName"] = magnetic.get_reference();
    device["Angle"] = 0;
    device["Disabled"] = false;
    device["Name"] = magnetic.get_reference();
    device["Id"] = generate_id();
    device["Parameters"]["Name"] = magnetic.get_reference();
    device["SubcircuitDefinition"]["Devices"] = ordered_json::array();
    device["SubcircuitDefinition"]["Connectors"] = ordered_json::array();
    device["SubcircuitDefinition"]["Name"] = magnetic.get_reference();
    device["SubcircuitDefinition"]["Id"] = generate_id();
    device["SubcircuitDefinition"]["Variables"] = ordered_json::array();
    device["SubcircuitDefinition"]["VariableFile"] = "";
    device["SubcircuitDefinitionID"] = device["SubcircuitDefinition"]["Id"];

    auto columns = core.get_columns();
    auto coreEffectiveArea = core.get_effective_area();

    auto coreEffectiveLengthMinusColumns = core.get_effective_length();
    if (columns.size() > 1) {
        for (auto column : columns) {
            if (column.get_coordinates()[0] >= 0) {
                coreEffectiveLengthMinusColumns -= column.get_height();
            }
        }
    }

    std::vector<std::vector<int>> columnBottomCoordinates; 
    std::vector<std::vector<int>> columnTopCoordinates;
    
    auto acResistanceCoefficientsPerWinding = CircuitSimulatorExporter::calculate_ac_resistance_coefficients_per_winding(magnetic, temperature);
    // Resolve AUTO mode from settings (fracpole support)
    auto resolvedMode_sb = CircuitSimulatorExporter::resolve_curve_fitting_mode(mode);
    std::vector<FractionalPoleNetwork> fracpoleNets_sb;
    if (resolvedMode_sb == CircuitSimulatorExporterCurveFittingModes::FRACPOLE) {
        fracpoleNets_sb = CircuitSimulatorExporter::calculate_fracpole_networks_per_winding(magnetic, temperature);
    }
    // Core losses using R-L ladder (calculate coefficients for positioning)
    auto coreResistanceCoefficients_sb = CircuitSimulatorExporter::calculate_core_resistance_coefficients(magnetic, temperature);
    double leakageInductance = resolve_dimensional_values(LeakageInductance().calculate_leakage_inductance(magnetic, frequency).get_leakage_inductance_per_winding()[0]);
    int numberLadderPairElements = acResistanceCoefficientsPerWinding[0].size() / 2 - 1;
    int numberCoreLadderPairElements = static_cast<int>(coreResistanceCoefficients_sb.size() / 2) + 1;

    for (size_t columnIndex = 0; columnIndex < columns.size(); ++columnIndex) {
        auto column = columns[columnIndex];
        auto gapsInThisColumn = core.find_gaps_by_column(column);
        ordered_json coreChunkJson;
        std::vector<int> columnCoordinates;
        if (column.get_coordinates()[0] == 0 && column.get_coordinates()[2] != 0) {
            columnCoordinates = {static_cast<int>(column.get_coordinates()[2] * _scale), 0};
        }
        else {
            columnCoordinates = {static_cast<int>(column.get_coordinates()[0] * _scale), 0};
        }

        std::vector<int> columnTopCoordinate = {columnCoordinates[0] + 3, -2 - numberLadderPairElements * 6};  // Don't ask
        std::vector<int> columnBottomCoordinate = {columnCoordinates[0] + 3, 4 + numberCoreLadderPairElements * 6};  // Don't ask
        int currentColumnGapHeight = 6;
        for (size_t gapIndex = 0; gapIndex < gapsInThisColumn.size(); ++gapIndex) {
            currentColumnGapHeight += 6;
            columnBottomCoordinate[1] += 6;
        }
        columnBottomCoordinates.push_back(columnBottomCoordinate);
        columnTopCoordinates.push_back(columnTopCoordinate);
    }


    std::vector<int> windingCoordinates = {columnTopCoordinates[0][0] - 2, columnTopCoordinates[0][1] - 6 + numberLadderPairElements * 6};
    windingCoordinates[1] -= 6 * (coil.get_functional_description().size() - 1);
    int numberLeftWindings = 0;
    int numberRightWindings = 0;
    int maximumWindingCoordinate = windingCoordinates[1];
    int maximumLadderCoordinate = windingCoordinates[1];
    int maximumLeftCoordinate = windingCoordinates[0];
    int maximumRightCoordinate = windingCoordinates[0];

    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        auto dcResistanceThisWinding = WindingLosses::calculate_effective_resistance_of_winding(magnetic, windingIndex, 0.1, temperature);
        auto winding = coil.get_functional_description()[windingIndex];
        std::vector<int> coordinates = windingCoordinates;
        ordered_json windingJson;
        ordered_json topPinJson;
        ordered_json bottomPinJson;
        ordered_json acResistorJson;
        std::vector<ordered_json> ladderJsons;
        std::vector<ordered_json> ladderConnectorsJsons;

        if (winding.get_isolation_side() == IsolationSide::PRIMARY) {
            windingJson = create_winding(winding.get_number_turns(), coordinates, 0, winding.get_name());
            coordinates[0] -= numberLeftWindings * 18;
        }
        else {
            windingJson = create_winding(winding.get_number_turns(), coordinates, 180, winding.get_name());
            coordinates[0] += numberRightWindings * 12;
        }

        if (winding.get_isolation_side() == IsolationSide::PRIMARY) {
            if (windingIndex == 0) {
                coordinates[0] -= 6;
                auto leakageInductanceJson = create_inductor(leakageInductance, coordinates, 0, winding.get_name() + " Leakage inductance");
                device["SubcircuitDefinition"]["Devices"].push_back(leakageInductanceJson);
            }
            coordinates[0] -= 6;
            acResistorJson = create_resistor(dcResistanceThisWinding, coordinates, 180, winding.get_name() + " AC resistance");

            std::pair<std::vector<ordered_json>, std::vector<ordered_json>> aux;
            if (resolvedMode_sb == CircuitSimulatorExporterCurveFittingModes::FRACPOLE &&
                windingIndex < fracpoleNets_sb.size()) {
                aux = create_fracpole_ladder(fracpoleNets_sb[windingIndex], coordinates, winding.get_name());
            } else {
                aux = create_ladder(acResistanceCoefficientsPerWinding[windingIndex], coordinates, winding.get_name());
            }
            if (acResistanceCoefficientsPerWinding[windingIndex].size() > 0) {
                coordinates[0] -= 6;
            }
            ladderJsons = aux.first;
            ladderConnectorsJsons = aux.second;

            coordinates[0] -= 2;
            bottomPinJson = create_pin(coordinates, 0, winding.get_name() + " Input");

            if (numberLeftWindings > 0) {
                auto connectorJson = create_connector({windingCoordinates[0] - numberLeftWindings * 18, windingCoordinates[1] + 1}, {windingCoordinates[0], windingCoordinates[1] + 1}, "Connector between DC resistance in winding " + std::to_string(windingIndex));
                device["SubcircuitDefinition"]["Connectors"].push_back(connectorJson);
            }
            topPinJson = create_pin({windingCoordinates[0] - 2, windingCoordinates[1] + 4}, 0, winding.get_name() + " Output");
        }
        else {
            coordinates[0] += 4;
            acResistorJson = create_resistor(dcResistanceThisWinding, coordinates, 180, winding.get_name() + " AC resistance");

            std::pair<std::vector<ordered_json>, std::vector<ordered_json>> aux;
            if (resolvedMode_sb == CircuitSimulatorExporterCurveFittingModes::FRACPOLE &&
                windingIndex < fracpoleNets_sb.size()) {
                aux = create_fracpole_ladder(fracpoleNets_sb[windingIndex], {coordinates[0] + 12, coordinates[1]}, winding.get_name());
            } else {
                aux = create_ladder(acResistanceCoefficientsPerWinding[windingIndex], {coordinates[0] + 12, coordinates[1]}, winding.get_name());
            }
            ladderJsons = aux.first;
            ladderConnectorsJsons = aux.second;
            if (acResistanceCoefficientsPerWinding[windingIndex].size() > 0) {
                coordinates[0] += 6;
            }

            coordinates[0] += 6;
            topPinJson = create_pin(coordinates, 180, winding.get_name() + " Input");

            if (numberRightWindings > 0) {
                auto connectorJson = create_connector({windingCoordinates[0] + numberRightWindings * 12 + 4, windingCoordinates[1] + 1}, {windingCoordinates[0] + 4, windingCoordinates[1] + 1}, "Connector between DC resistance in winding " + std::to_string(windingIndex));
                device["SubcircuitDefinition"]["Connectors"].push_back(connectorJson);
            }
            bottomPinJson = create_pin({windingCoordinates[0] + 4, windingCoordinates[1] + 4}, 180, winding.get_name() + " Output");
        }

        std::vector<int> connectorTopCoordinates = {windingCoordinates[0] + 2, windingCoordinates[1] + 5};
        std::vector<int> connectorBottomCoordinates = {windingCoordinates[0] + 2, windingCoordinates[1] + 7};
        if (windingIndex == coil.get_functional_description().size() - 1) {
            connectorBottomCoordinates[1] = windingCoordinates[1] + 6;
        }

        auto connectorJson = create_connector(connectorBottomCoordinates, connectorTopCoordinates, "Connector between winding " + std::to_string(windingIndex) + " and winding " + std::to_string(windingIndex + 1));
        device["SubcircuitDefinition"]["Connectors"].push_back(connectorJson);

        device["SubcircuitDefinition"]["Devices"].push_back(windingJson);
        device["SubcircuitDefinition"]["Devices"].push_back(topPinJson);
        device["SubcircuitDefinition"]["Devices"].push_back(bottomPinJson);
        device["SubcircuitDefinition"]["Devices"].push_back(acResistorJson);
        for (auto ladderJson : ladderJsons) {
            device["SubcircuitDefinition"]["Devices"].push_back(ladderJson);
        }
        for (auto ladderConnectorsJson : ladderConnectorsJsons) {
            device["SubcircuitDefinition"]["Connectors"].push_back(ladderConnectorsJson);
        }

        windingCoordinates[1] += 6;
        if (winding.get_isolation_side() == IsolationSide::PRIMARY) {
            numberLeftWindings++;
        }
        else {
            numberRightWindings++;
        }
    }


    // Magnetizing current and core losses using R-L ladder
    // The R-L ladder provides impedance that INCREASES with frequency, matching core loss behavior.
    {
        auto coreResistanceCoefficients = CircuitSimulatorExporter::calculate_core_resistance_coefficients(magnetic, temperature);
        int numberCoreLadderPairElements = static_cast<int>(coreResistanceCoefficients.size() / 2);
        
        std::vector<int> coordinates = windingCoordinates;
        std::vector<ordered_json> ladderJsons;
        std::vector<ordered_json> ladderConnectorsJsons;
        coordinates[1] += numberCoreLadderPairElements * 6 - 5;

        {
            std::vector<int> connectorTopCoordinates = {windingCoordinates[0] + 2, windingCoordinates[1]};
            std::vector<int> connectorBottomCoordinates = {coordinates[0] + 2, coordinates[1] + 5};

            auto connectorJson = create_connector(connectorTopCoordinates, connectorBottomCoordinates, "Central column connector to core losses");
            device["SubcircuitDefinition"]["Connectors"].push_back(connectorJson);
        }

        // Add core loss R-L ladder using create_ladder
        if (!coreResistanceCoefficients.empty()) {
            coordinates[0] -= 6;
            auto aux = create_ladder(coreResistanceCoefficients, coordinates, "Core losses");
            if (coreResistanceCoefficients.size() > 0) {
                coordinates[0] -= 6;
            }
            ladderJsons = aux.first;
            ladderConnectorsJsons = aux.second;

            for (auto ladderJson : ladderJsons) {
                device["SubcircuitDefinition"]["Devices"].push_back(ladderJson);
            }
            for (auto ladderConnectorsJson : ladderConnectorsJsons) {
                device["SubcircuitDefinition"]["Connectors"].push_back(ladderConnectorsJson);
            }
        }
    }

    for (auto deviceJson : device["SubcircuitDefinition"]["Devices"]) {
        maximumLadderCoordinate = std::min(maximumLadderCoordinate, deviceJson["Top"].get<int>());
        maximumLeftCoordinate = std::max(maximumLeftCoordinate, deviceJson["Left"].get<int>());
        maximumRightCoordinate = std::min(maximumRightCoordinate, deviceJson["Left"].get<int>());
    }

    for (size_t columnIndex = 0; columnIndex < columns.size(); ++columnIndex) {
        columnTopCoordinates[columnIndex][1] = maximumLadderCoordinate - 5;
    }
    if (columns.size() > 1) {
        columnTopCoordinates[1][0] = maximumLeftCoordinate + 6;
        columnBottomCoordinates[1][0] = maximumLeftCoordinate + 6;
    }
    if (columns.size() > 2) {
        columnTopCoordinates[2][0] = maximumRightCoordinate - 2;
        columnBottomCoordinates[2][0] = maximumRightCoordinate - 2;
    }


    for (size_t columnIndex = 0; columnIndex < columns.size(); ++columnIndex) {
        auto column = columns[columnIndex];
        auto gapsInThisColumn = core.find_gaps_by_column(column);
        ordered_json coreChunkJson;
        std::vector<int> coordinates = columnBottomCoordinates[columnIndex];
        coordinates[1] -= (gapsInThisColumn.size() + 1) * 6 - 2;
        coordinates[0] -= 3;
        if (columnIndex == 0) {
            coreChunkJson = create_core(core.get_initial_permeability(), coordinates, coreEffectiveArea, coreEffectiveLengthMinusColumns, 90, "Core winding column and plates " + std::to_string(columnIndex));
        }
        else {
            coreChunkJson = create_core(core.get_initial_permeability(), coordinates, coreEffectiveArea, column.get_height(), 90, "Core lateral column " + std::to_string(columnIndex));

            std::vector<int> connectorTopCoordinates = {columnTopCoordinates[0][0], columnTopCoordinates[columnIndex][1]};
            std::vector<int> connectorBottomCoordinates = {columnTopCoordinates[columnIndex][0], coordinates[1] - 2};
            auto connectorJson = create_connector(connectorTopCoordinates, connectorBottomCoordinates, "Connector between column " + std::to_string(columnIndex) + " and top");
            device["SubcircuitDefinition"]["Connectors"].push_back(connectorJson);
        }
        device["SubcircuitDefinition"]["Devices"].push_back(coreChunkJson);

        for (size_t gapIndex = 0; gapIndex < gapsInThisColumn.size(); ++gapIndex) {
            coordinates[1] += 6;
            auto gap = gapsInThisColumn[gapIndex];
            if (!gap.get_coordinates()) {
                throw GapException("Gap is not processed");
            }
            std::vector<int> gapCoordinates = {coordinates[0], coordinates[1]};

            if (gap.get_length() > 0) {

                auto gapJson = create_air_gap(gapCoordinates, gap.get_area().value(), gap.get_length(), 90, "Column " + std::to_string(columnIndex) + " gap " + std::to_string(gapIndex));
                device["SubcircuitDefinition"]["Devices"].push_back(gapJson);
            }
            else {
                std::vector<int> zeroGapConnectorTopCoordinates = {gapCoordinates[0] + 3, gapCoordinates[1]- 2};
                std::vector<int> zeroGapConnectorBottomCoordinates = {gapCoordinates[0] + 3, gapCoordinates[1] + 4};
                auto connectorJson = create_connector(zeroGapConnectorTopCoordinates, zeroGapConnectorBottomCoordinates, "Connector replacing gap of 0 length");
                device["SubcircuitDefinition"]["Connectors"].push_back(connectorJson);
            }

        }
    }

    {
        std::vector<int> finalConnectorTopCoordinates = {columnTopCoordinates[0][0], columnTopCoordinates[0][1]};
        std::vector<int> finalConnectorBottomCoordinates = {columnTopCoordinates[0][0], maximumWindingCoordinate + 1};
        auto connectorJson = create_connector(finalConnectorTopCoordinates, finalConnectorBottomCoordinates, "Connector between winding 0 and top");
        device["SubcircuitDefinition"]["Connectors"].push_back(connectorJson);
    }

    // Magnetic ground
    {
        std::vector<int> columnBottomCoordinatesAux = {0, columnTopCoordinates[0][1]};
        columnBottomCoordinatesAux[0] += 2;
        columnBottomCoordinatesAux[1] -= 2;
        auto groundJson = create_magnetic_ground(columnBottomCoordinatesAux, 180, "Magnetic ground");
        device["SubcircuitDefinition"]["Devices"].push_back(groundJson);
    }

    // Horizontal bottom connectors
    if (columns.size() == 1) {
        auto columnBottomCoordinatesAux = columnBottomCoordinates[0];
        columnBottomCoordinatesAux[1] = 0;
        columnBottomCoordinatesAux[0] += _modelSize / 2;
        auto bottomConnectorJson = create_connector(columnBottomCoordinates[0], columnBottomCoordinatesAux, "Bottom Connector between column " + std::to_string(0) + " and middle");
        device["SubcircuitDefinition"]["Connectors"].push_back(bottomConnectorJson);
        auto topConnectorJson = create_connector(columnTopCoordinates[0], columnBottomCoordinatesAux, "Top Connector between column " + std::to_string(0) + " and middle");
        device["SubcircuitDefinition"]["Connectors"].push_back(topConnectorJson);
        device["SubcircuitDefinition"]["Connectors"] = merge_connectors(device["SubcircuitDefinition"]["Connectors"]);
    }
    else {
        for (size_t columnIndex = 1; columnIndex < columns.size(); ++columnIndex) {
            auto connectorJson = create_connector(columnBottomCoordinates[0], columnBottomCoordinates[columnIndex], "Bottom Connector between column " + std::to_string(0) + " and columm " + std::to_string(columnIndex));
            device["SubcircuitDefinition"]["Connectors"].push_back(connectorJson);
        }
    }
    device["SubcircuitDefinition"]["Connectors"] = merge_connectors(device["SubcircuitDefinition"]["Connectors"]);

    library["Devices"] = ordered_json::array();
    library["Devices"].push_back(device);
    simulation["Libraries"].push_back(library);
    return simulation.dump(2);
}

// =============================================================================
// FRACPOLE SPICE EMITTERS
// =============================================================================

// Emit the fracpole subcircuit + instantiation for winding i in NgSpice/LTSpice format.
// The winding resistance model is:
//   Rdc in series with the skin_effect subcircuit (fracpole + gyrator).
//   P+ ---[Rdc]---[skin_effect]--- Node_R_Lmag
static std::string emit_fracpole_winding_spice(
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
// Coefficients are [R1, L1, R2, L2, R3, L3] from calculate_core_resistance_coefficients().
static std::string emit_core_ladder_spice(
    const std::vector<double>& coeffs, size_t numWindings) {
    
    if (coeffs.size() < 2) {
        return "";  // Not enough coefficients
    }
    
    std::string s;
    s += "* Core loss R-L ladder (" + std::to_string(coeffs.size()/2) + " stages)\n";
    
    // Build ladder: series of (L || R) stages
    // Z = L1 || (R1 + L2 || (R2 + L3 || R3))
    // Start from innermost stage and work outward
    std::string lastNode = "P1-";
    int stageCount = static_cast<int>(coeffs.size() / 2);
    
    for (int k = stageCount - 1; k >= 0; --k) {
        double R = coeffs[k * 2];
        double L = coeffs[k * 2 + 1];
        std::string kk = std::to_string(k + 1);
        std::string thisNode = (k == 0) ? "Node_R_Lmag_1" : ("Node_Rcore_" + std::to_string(k));
        
        // Validate values
        if (R <= 0 || L <= 0 || R > 1e6 || L > 1) {
            continue;  // Skip invalid stage
        }
        
        // R in series path
        s += "Rcore" + kk + " " + thisNode + " Node_Lcore_" + kk + " " + to_string(R, 12) + "\n";
        // L in parallel with the rest (connects to lastNode)
        s += "Lcore" + kk + " " + thisNode + " " + lastNode + " " + to_string(L, 12) + "\n";
        
        lastNode = "Node_Lcore_" + kk;
    }
    
    // Connect the final node to ground (P1-)
    if (lastNode != "P1-") {
        s += "Rcore_final " + lastNode + " P1- 1e-6\n";  // Near-zero resistance to close the circuit
    }
    
    return s;
}

// Emit saturating magnetizing inductance for ngspice using behavioral modeling.
// Uses a tanh-based flux linkage model: Lambda(I) = Lambda_sat * tanh(I/I_sat) + L_linear * I
// where L_linear accounts for the air gap contribution.
static std::string emit_saturating_inductor_ngspice(
    const SaturationParameters& sat,
    const std::string& windingIndex,
    const std::string& nodeIn,
    const std::string& nodeOut) {
    
    if (!sat.valid) {
        return "";  // Fall back to linear inductor
    }
    
    std::string s;
    const double mu0 = 4e-7 * M_PI;
    
    // Calculate saturation current and flux linkage
    double Isat = sat.Isat();
    double lambdaSat = sat.fluxLinkageSat();
    
    // The model uses: V = d(Lambda)/dt where Lambda = Lambda_sat * tanh(I/I_sat) + L_gap * I
    // L_gap is the linear component (from air gap), approximately L_mag * (1 - 1/effective_permeability)
    // For a typical gapped core, most inductance comes from the gap, so L_gap ≈ L_mag
    double effectiveMuR = sat.Lmag / (mu0 * sat.primaryTurns * sat.primaryTurns * sat.Ae / sat.le);
    double gapFactor = 1.0 / effectiveMuR;  // Fraction of reluctance from gap
    double Lgap = sat.Lmag * gapFactor;  // Linear (gap) component
    
    // Ensure reasonable values
    if (Isat < 1e-6) Isat = 1e-6;
    if (lambdaSat < 1e-9) lambdaSat = 1e-9;
    if (Lgap < 1e-12) Lgap = sat.Lmag * 0.5;  // Fallback
    
    s += "* Saturating magnetizing inductance (winding " + windingIndex + ")\n";
    s += "* Bsat=" + to_string(sat.Bsat, 4) + "T, Isat=" + to_string(Isat, 4) + "A\n";
    
    // ngspice behavioral inductor using flux linkage
    // V = dLambda/dt = dLambda/dI * dI/dt
    // For Lambda = Lsat*tanh(I/Isat) + Lgap*I:
    // dLambda/dI = Lsat/Isat * sech²(I/Isat) + Lgap
    //            = Lsat/Isat * (1 - tanh²(I/Isat)) + Lgap
    // We use the GVALUE element with Laplace to model this
    
    // Simpler approach: Use the inductor with polynomial current dependence
    // L(I) = L0 / (1 + (I/Isat)^2) gives smooth saturation
    double L0 = sat.Lmag;
    
    s += ".param Lmag_" + windingIndex + "_L0=" + to_string(L0, 12) + "\n";
    s += ".param Lmag_" + windingIndex + "_Isat=" + to_string(Isat, 6) + "\n";
    
    // Use behavioral source: V = L(I) * dI/dt
    // L(I) = L0 / (1 + (abs(I)/Isat)^2)
    // Note: ngspice behavioral sources use I() for current through voltage source
    std::string senseName = "Vmag_sense_" + windingIndex;
    std::string fluxName = "Lmag_flux_" + windingIndex;
    
    // Current sense (zero volt source)
    s += senseName + " " + nodeIn + " " + fluxName + " 0\n";
    
    // Behavioral inductor: V = L(I) * dI/dt
    // Using EL (behavioral voltage source) with ddt() function
    s += "BLmag_" + windingIndex + " " + fluxName + " " + nodeOut;
    s += " V='Lmag_" + windingIndex + "_L0/(1+pow(abs(I(" + senseName + "))/Lmag_" + windingIndex + "_Isat,2))*ddt(I(" + senseName + "))'\n";
    
    return s;
}

// Emit saturating magnetizing inductance for LTspice using Flux= syntax.
// LTspice supports nonlinear inductors with Flux={expression} directly.
static std::string emit_saturating_inductor_ltspice(
    const SaturationParameters& sat,
    const std::string& windingIndex,
    const std::string& nodeIn,
    const std::string& nodeOut) {
    
    if (!sat.valid) {
        return "";  // Fall back to linear inductor
    }
    
    std::string s;
    const double mu0 = 4e-7 * M_PI;
    
    // Calculate saturation current and flux linkage
    double Isat = sat.Isat();
    double lambdaSat = sat.fluxLinkageSat();
    
    // Calculate gap contribution
    double effectiveMuR = sat.Lmag / (mu0 * sat.primaryTurns * sat.primaryTurns * sat.Ae / sat.le);
    double gapFactor = std::min(0.99, 1.0 / effectiveMuR);
    double Lgap = sat.Lmag * gapFactor;
    
    // Ensure reasonable values  
    if (Isat < 1e-6) Isat = 1e-6;
    if (lambdaSat < 1e-9) lambdaSat = 1e-9;
    
    s += "* Saturating magnetizing inductance (winding " + windingIndex + ")\n";
    s += "* Bsat=" + to_string(sat.Bsat, 4) + "T, Isat=" + to_string(Isat, 4) + "A\n";
    
    // LTspice flux linkage model: Lambda = Lambda_sat * tanh(I/I_sat) + L_gap * I
    // The tanh term models core saturation, L_gap*I models the linear air gap
    s += ".param Lambda_sat_" + windingIndex + "=" + to_string(lambdaSat, 12) + "\n";
    s += ".param I_sat_" + windingIndex + "=" + to_string(Isat, 6) + "\n";
    s += ".param L_gap_" + windingIndex + "=" + to_string(Lgap, 12) + "\n";
    
    // Use LTspice Flux= syntax for nonlinear inductor
    // x is the current through the inductor
    s += "Lmag_" + windingIndex + " " + nodeIn + " " + nodeOut;
    s += " Flux={Lambda_sat_" + windingIndex + "*tanh(x/I_sat_" + windingIndex + ")+L_gap_" + windingIndex + "*x}\n";
    
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

static std::string emit_mutual_resistance_network_spice(
    const std::vector<CircuitSimulatorExporter::MutualResistanceCoefficients>& mutualCoeffs,
    double magnetizingInductance,
    size_t numWindings) {
    
    if (mutualCoeffs.empty() || numWindings < 2) {
        return "";  // No mutual resistance to model
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

std::string CircuitSimulatorExporterNgspiceModel::export_magnetic_as_subcircuit(Magnetic magnetic, double frequency, double temperature, std::optional<std::string> filePathOrFile, CircuitSimulatorExporterCurveFittingModes mode) {
    std::string headerString = "* Magnetic model made with OpenMagnetics\n";
    headerString += "* " + magnetic.get_reference() + "\n\n";
    headerString += ".subckt " + fix_filename(magnetic.get_reference());
    std::string circuitString = "";
    std::string parametersString = "";
    std::string footerString = ".ends " + fix_filename(magnetic.get_reference());

    auto coil = magnetic.get_coil();

    double magnetizingInductance = resolve_dimensional_values(MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(magnetic).get_magnetizing_inductance());
    auto acResistanceCoefficientsPerWinding = CircuitSimulatorExporter::calculate_ac_resistance_coefficients_per_winding(magnetic, temperature);
    auto leakageInductances = LeakageInductance().calculate_leakage_inductance(magnetic, Defaults().measurementFrequency).get_leakage_inductance_per_winding();

    // Resolve AUTO mode from settings
    auto resolvedMode = CircuitSimulatorExporter::resolve_curve_fitting_mode(mode);
    std::vector<FractionalPoleNetwork> fracpoleNets;
    std::optional<FractionalPoleNetwork> coreFracNet;
    if (resolvedMode == CircuitSimulatorExporterCurveFittingModes::FRACPOLE) {
        fracpoleNets = CircuitSimulatorExporter::calculate_fracpole_networks_per_winding(magnetic, temperature);
        try { coreFracNet = CircuitSimulatorExporter::calculate_core_fracpole_network(magnetic, temperature); }
        catch (...) {}
    }

    parametersString += ".param MagnetizingInductance_Value=" + std::to_string(magnetizingInductance) + "\n";
    parametersString += ".param Permeance=MagnetizingInductance_Value/NumberTurns_1**2\n";
    
    // Check if saturation modeling is enabled
    auto& settings = Settings::GetInstance();
    bool includeSaturation = settings.get_circuit_simulator_include_saturation();
    SaturationParameters satParams;
    if (includeSaturation) {
        satParams = get_saturation_parameters(magnetic, temperature);
        // Saturation only supported for single winding (inductor) - K coupling doesn't work with behavioral sources
        if (satParams.valid && coil.get_functional_description().size() > 1) {
            // For transformers, fall back to linear model with warning comment
            circuitString += "* WARNING: Saturation modeling not supported for multi-winding transformers in ngspice\n";
            circuitString += "* K coupling requires linear inductors. Using linear model.\n";
            includeSaturation = false;
        }
    }
    
    // Store coupling coefficients for each pair (indexed from winding 1)
    std::vector<double> couplingCoeffs;
    size_t numWindings = coil.get_functional_description().size();
    
    for (size_t index = 0; index < numWindings; index++) {
        auto effectiveResistanceThisWinding = WindingLosses::calculate_effective_resistance_of_winding(magnetic, index, 0.1, temperature);
        std::string is = std::to_string(index + 1);
        parametersString += ".param Rdc_" + is + "_Value=" + std::to_string(effectiveResistanceThisWinding) + "\n";
        parametersString += ".param NumberTurns_" + is + "=" + std::to_string(coil.get_functional_description()[index].get_number_turns()) + "\n";
        if (index > 0) {
            double leakageInductance = resolve_dimensional_values(leakageInductances[index - 1]);
            if (leakageInductance >= magnetizingInductance) {
                leakageInductance = magnetizingInductance * 0.1;
            }
            double couplingCoefficient = sqrt((magnetizingInductance - leakageInductance) / magnetizingInductance);
            couplingCoeffs.push_back(couplingCoefficient);
            parametersString += ".param Llk_" + is + "_Value=" + std::to_string(leakageInductance) + "\n";
            parametersString += ".param CouplingCoefficient_1" + is + "_Value=" + std::to_string(couplingCoefficient) + "\n";
        }

        std::vector<std::string> c = to_string(acResistanceCoefficientsPerWinding[index]);

        // FRACPOLE mode: emit fracpole-based winding subcircuit
        if (resolvedMode == CircuitSimulatorExporterCurveFittingModes::FRACPOLE) {
            if (index < fracpoleNets.size()) {
                circuitString += emit_fracpole_winding_spice(fracpoleNets[index], index, is, effectiveResistanceThisWinding);
            } else {
                circuitString += "Rdc" + is + " P" + is + "+ Node_R_Lmag_" + is + " {Rdc_" + is + "_Value}\n";
            }
            // Emit magnetizing inductance (saturating or linear)
            if (includeSaturation && satParams.valid) {
                circuitString += emit_saturating_inductor_ngspice(satParams, is, "Node_R_Lmag_" + is, "P" + is + "-");
            } else {
                circuitString += "Lmag_" + is + " Node_R_Lmag_" + is + " P" + is + "- {NumberTurns_" + is + "**2*Permeance}\n";
            }
        }
        else if (resolvedMode == CircuitSimulatorExporterCurveFittingModes::ANALYTICAL) {
            throw std::invalid_argument("Analytical mode not supported in NgSpice");
        }
        else {
            // LADDER mode (default)
            bool validLadderCoeffs = true;
            for (size_t ladderIndex = 0; ladderIndex < acResistanceCoefficientsPerWinding[index].size(); ladderIndex+=2) {
                double inductanceVal = acResistanceCoefficientsPerWinding[index][ladderIndex + 1];
                double resistanceVal = acResistanceCoefficientsPerWinding[index][ladderIndex];
                if (inductanceVal > 0.1 || resistanceVal > 100.0 || inductanceVal < 0 || resistanceVal < 0) {
                    validLadderCoeffs = false;
                    break;
                }
            }
            
            if (validLadderCoeffs && acResistanceCoefficientsPerWinding[index].size() >= 2) {
                for (size_t ladderIndex = 0; ladderIndex < acResistanceCoefficientsPerWinding[index].size(); ladderIndex+=2) {
                    std::string ladderIndexs = std::to_string(ladderIndex);
                    circuitString += "Lladder" + is + "_" + ladderIndexs + " P" + is + "+ Node_Lladder_" + is + "_" + ladderIndexs + " " + c[ladderIndex + 1] + "\n";
                    if (ladderIndex == 0) {
                        circuitString += "Rladder" + is + "_" + ladderIndexs + " Node_Lladder_" + is + "_" + ladderIndexs + " Node_R_Lmag_" + is + " " + c[ladderIndex] + "\n";
                    }
                    else {
                        circuitString += "Rladder" + is + "_" + ladderIndexs + " Node_Lladder_" + is + "_" + ladderIndexs + " Node_Lladder_" + is + "_" + std::to_string(ladderIndex - 2) + " " + c[ladderIndex] + "\n";
                    }
                }
                circuitString += "Rdc" + is + " P" + is + "+ Node_R_Lmag_" + is + " {Rdc_" + is + "_Value}\n";
            } else {
                circuitString += "Rdc" + is + " P" + is + "+ Node_R_Lmag_" + is + " {Rdc_" + is + "_Value}\n";
            }
            // Emit magnetizing inductance (saturating or linear)
            if (includeSaturation && satParams.valid) {
                circuitString += emit_saturating_inductor_ngspice(satParams, is, "Node_R_Lmag_" + is, "P" + is + "-");
            } else {
                circuitString += "Lmag_" + is + " Node_R_Lmag_" + is + " P" + is + "- {NumberTurns_" + is + "**2*Permeance}\n";
            }
        }

        headerString += " P" + is + "+ P" + is + "-";
    }
    
    // Generate pairwise K statements for magnetic coupling
    // ngspice has a bug with 3+ inductor K statements inside subcircuits - the third inductor
    // doesn't get the proper hierarchical path prefix. Use pairwise K statements as workaround.
    // Each K statement gets a unique name (K12, K13, K23, etc.)
    // Use per-pair leakage inductance calculation for accurate coupling coefficients
    if (numWindings == 2) {
        // Simple 2-winding case - use already calculated coupling, capped at 0.98 for stability
        double k12 = couplingCoeffs.size() > 0 ? std::min(0.98, couplingCoeffs[0]) : 0.98;
        circuitString += "K Lmag_1 Lmag_2 " + std::to_string(k12) + "\n";
    } else if (numWindings >= 3) {
        // For 3+ windings, calculate coupling for each pair
        for (size_t i = 0; i < numWindings; i++) {
            for (size_t j = i + 1; j < numWindings; j++) {
                // Calculate leakage inductance between winding i and j
                double leakageIJ = 0.0;
                try {
                    auto leakageResult = LeakageInductance().calculate_leakage_inductance(
                        magnetic, Defaults().measurementFrequency, i, j);
                    auto leakagePerWinding = leakageResult.get_leakage_inductance_per_winding();
                    if (leakagePerWinding.size() > 0) {
                        leakageIJ = resolve_dimensional_values(leakagePerWinding[0]);
                    }
                } catch (...) {
                    // Fallback: use high coupling if calculation fails
                    leakageIJ = magnetizingInductance * 0.02;  // ~99% coupling
                }
                
                // Clamp leakage inductance to avoid negative or very low coupling
                if (leakageIJ >= magnetizingInductance) {
                    leakageIJ = magnetizingInductance * 0.1;  // Limit to ~95% coupling minimum
                }
                if (leakageIJ < 0) {
                    leakageIJ = magnetizingInductance * 0.02;  // ~99% coupling
                }
                
                double kij = sqrt((magnetizingInductance - leakageIJ) / magnetizingInductance);
                // Ensure coupling is in valid range
                // Cap at 0.98 for numerical stability (matches IDEAL mode default)
                kij = std::max(0.5, std::min(0.98, kij));
                
                std::string kName = "K" + std::to_string(i + 1) + std::to_string(j + 1);
                circuitString += kName + " Lmag_" + std::to_string(i + 1) + " Lmag_" + std::to_string(j + 1) + " " + std::to_string(kij) + "\n";
            }
        }
    }

    // Core losses: Use R-L ladder network
    // The R-L ladder provides impedance that INCREASES with frequency, which matches
    // the physical behavior of core losses (higher losses at higher frequencies).
    // This is opposite to the R-C fracpole network which decreases with frequency.
    auto coreResistanceCoefficients = CircuitSimulatorExporter::calculate_core_resistance_coefficients(magnetic, temperature);
    if (!coreResistanceCoefficients.empty()) {
        circuitString += emit_core_ladder_spice(coreResistanceCoefficients, numWindings);
    }

    // Mutual resistance: auxiliary winding network for cross-coupling losses (Hesterman 2020)
    // Only for multi-winding transformers when setting is enabled
    bool includeMutualResistance = settings.get_circuit_simulator_include_mutual_resistance();
    if (includeMutualResistance && numWindings >= 2) {
        auto mutualResistanceCoeffs = CircuitSimulatorExporter::calculate_mutual_resistance_coefficients(magnetic, temperature);
        if (!mutualResistanceCoeffs.empty()) {
            circuitString += emit_mutual_resistance_network_spice(mutualResistanceCoeffs, magnetizingInductance, numWindings);
        }
    }

    return headerString + "\n" + circuitString + "\n" + parametersString + "\n" + footerString;
}

std::string CircuitSimulatorExporterLtspiceModel::export_magnetic_as_subcircuit(Magnetic magnetic, double frequency, double temperature, std::optional<std::string> filePathOrFile, CircuitSimulatorExporterCurveFittingModes mode) {
    std::string headerString = "* Magnetic model made with OpenMagnetics\n";
    headerString += "* " + magnetic.get_reference() + "\n\n";
    headerString += ".subckt " + fix_filename(magnetic.get_reference());
    std::string circuitString = "";
    std::string parametersString = "";
    std::string footerString = ".ends " + fix_filename(magnetic.get_reference());

    auto coil = magnetic.get_coil();

    double magnetizingInductance = resolve_dimensional_values(MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(magnetic).get_magnetizing_inductance());
    auto acResistanceCoefficientsPerWinding = CircuitSimulatorExporter::calculate_ac_resistance_coefficients_per_winding(magnetic, temperature);
    auto leakageInductances = LeakageInductance().calculate_leakage_inductance(magnetic, Defaults().measurementFrequency).get_leakage_inductance_per_winding();

    // Resolve AUTO mode from settings
    auto resolvedMode_lt = CircuitSimulatorExporter::resolve_curve_fitting_mode(mode);
    std::vector<FractionalPoleNetwork> fracpoleNets_lt;
    std::optional<FractionalPoleNetwork> coreFracNet_lt;
    if (resolvedMode_lt == CircuitSimulatorExporterCurveFittingModes::FRACPOLE) {
        fracpoleNets_lt = CircuitSimulatorExporter::calculate_fracpole_networks_per_winding(magnetic, temperature);
        try { coreFracNet_lt = CircuitSimulatorExporter::calculate_core_fracpole_network(magnetic, temperature); }
        catch (...) {}
    }

    parametersString += ".param MagnetizingInductance_Value=" + std::to_string(magnetizingInductance) + "\n";
    parametersString += ".param Permeance=MagnetizingInductance_Value/NumberTurns_1**2\n";
    
    // Check if saturation modeling is enabled
    auto& settings = Settings::GetInstance();
    bool includeSaturation = settings.get_circuit_simulator_include_saturation();
    SaturationParameters satParams;
    if (includeSaturation) {
        satParams = get_saturation_parameters(magnetic, temperature);
        // Saturation only supported for single winding (inductor) - K coupling doesn't work with nonlinear inductors
        if (satParams.valid && coil.get_functional_description().size() > 1) {
            // For transformers, fall back to linear model with warning comment
            circuitString += "* WARNING: Saturation modeling not supported for multi-winding transformers in LTspice\n";
            circuitString += "* K coupling requires linear inductors. Using linear model.\n";
            includeSaturation = false;
        }
    }
    
    for (size_t index = 0; index < coil.get_functional_description().size(); index++) {
        auto effectiveResistanceThisWinding = WindingLosses::calculate_effective_resistance_of_winding(magnetic, index, 0.1, temperature);
        std::string is = std::to_string(index + 1);
        parametersString += ".param Rdc_" + is + "_Value=" + std::to_string(effectiveResistanceThisWinding) + "\n";
        parametersString += ".param NumberTurns_" + is + "=" + std::to_string(coil.get_functional_description()[index].get_number_turns()) + "\n";
        if (index > 0) {
            double leakageInductance = resolve_dimensional_values(leakageInductances[index - 1]);
            if (leakageInductance >= magnetizingInductance) {
                leakageInductance = magnetizingInductance * 0.1;
            }
            double couplingCoefficient = sqrt((magnetizingInductance - leakageInductance) / magnetizingInductance);
            parametersString += ".param Llk_" + is + "_Value=" + std::to_string(leakageInductance) + "\n";
            parametersString += ".param CouplingCoefficient_1" + is + "_Value=" + std::to_string(couplingCoefficient) + "\n";
        }

        std::vector<std::string> c = to_string(acResistanceCoefficientsPerWinding[index]);

        // FRACPOLE mode: emit fracpole-based winding subcircuit
        if (resolvedMode_lt == CircuitSimulatorExporterCurveFittingModes::FRACPOLE) {
            if (index < fracpoleNets_lt.size()) {
                circuitString += emit_fracpole_winding_spice(fracpoleNets_lt[index], index, is, effectiveResistanceThisWinding);
            } else {
                circuitString += "Rdc" + is + " P" + is + "+ Node_R_Lmag_" + is + " {Rdc_" + is + "_Value}\n";
            }
            // Emit magnetizing inductance (saturating or linear)
            if (includeSaturation && satParams.valid) {
                circuitString += emit_saturating_inductor_ltspice(satParams, is, "Node_R_Lmag_" + is, "P" + is + "-");
            } else {
                circuitString += "Lmag_" + is + " Node_R_Lmag_" + is + " P" + is + "- {NumberTurns_" + is + "**2*Permeance}\n";
            }
        }
        else if (resolvedMode_lt == CircuitSimulatorExporterCurveFittingModes::ANALYTICAL) {
            circuitString += "E" + is + " P" + is + "+ Node_R_Lmag_" + is + " P" + is + "+ Node_R_Lmag_" + is + " Laplace = 1 /(" + c[0] + " + " + c[1] + " * sqrt(abs(s)/(2*pi)) + " + c[2] + " * abs(s)/(2*pi))\n";
            // Emit magnetizing inductance (saturating or linear)
            if (includeSaturation && satParams.valid) {
                circuitString += emit_saturating_inductor_ltspice(satParams, is, "P" + is + "-", "Node_R_Lmag_" + is);
            } else {
                circuitString += "Lmag_" + is + " P" + is + "- Node_R_Lmag_" + is + " {NumberTurns_" + is + "**2*Permeance}\n";
            }
        }
        else {
            // LADDER mode (default)
            // Check if ladder coefficients are valid (inductance values should be small, < 0.1H)
            // If fitting failed, coefficients can be very large (1.0H or more) which breaks the circuit
            bool validLadderCoeffs = true;
            for (size_t ladderIndex = 0; ladderIndex < acResistanceCoefficientsPerWinding[index].size(); ladderIndex+=2) {
                double inductanceVal = acResistanceCoefficientsPerWinding[index][ladderIndex + 1];
                double resistanceVal = acResistanceCoefficientsPerWinding[index][ladderIndex];
                // Sanity check: inductance should be < 0.1H (100mH), resistance should be < 100 Ohm per ladder element
                if (inductanceVal > 0.1 || resistanceVal > 100.0 || inductanceVal < 0 || resistanceVal < 0) {
                    validLadderCoeffs = false;
                    break;
                }
            }
            
            if (validLadderCoeffs && acResistanceCoefficientsPerWinding[index].size() >= 2) {
                // Use full ladder network for AC resistance modeling
                for (size_t ladderIndex = 0; ladderIndex < acResistanceCoefficientsPerWinding[index].size(); ladderIndex+=2) {
                    std::string ladderIndexs = std::to_string(ladderIndex);
                    circuitString += "Lladder" + is + "_" + ladderIndexs + " P" + is + "+ Node_Lladder_" + is + "_" + ladderIndexs + " " + c[ladderIndex + 1] + "\n";
                    if (ladderIndex == 0) {
                        circuitString += "Rladder" + is + "_" + ladderIndexs + " Node_Lladder_" + is + "_" + ladderIndexs + " Node_R_Lmag_" + is + " " + c[ladderIndex] + "\n";
                    }
                    else {
                        circuitString += "Rladder" + is + "_" + ladderIndexs + " Node_Lladder_" + is + "_" + ladderIndexs + " Node_Lladder_" + is + "_" + std::to_string(ladderIndex - 2) + " " + c[ladderIndex] + "\n";
                    }
                }
                circuitString += "Rdc" + is + " P" + is + "+ Node_R_Lmag_" + is + " {Rdc_" + is + "_Value}\n";
            } else {
                // Ladder fitting failed - use simplified model with just DC resistance
                circuitString += "Rdc" + is + " P" + is + "+ Node_R_Lmag_" + is + " {Rdc_" + is + "_Value}\n";
            }
            // Emit magnetizing inductance (saturating or linear)
            if (includeSaturation && satParams.valid) {
                circuitString += emit_saturating_inductor_ltspice(satParams, is, "Node_R_Lmag_" + is, "P" + is + "-");
            } else {
                circuitString += "Lmag_" + is + " Node_R_Lmag_" + is + " P" + is + "- {NumberTurns_" + is + "**2*Permeance}\n";
            }
        }
        if (index > 0) {
            // Each K statement needs a unique name in LTspice (K1, K2, etc.)
            circuitString += "K" + is + " Lmag_1 Lmag_" + is + " {CouplingCoefficient_1" + is + "_Value}\n";
        }

        headerString += " P" + is + "+ P" + is + "-";

    }

    // Core losses: Use R-L ladder network (same as ngspice)
    auto coreResistanceCoefficients = CircuitSimulatorExporter::calculate_core_resistance_coefficients(magnetic, temperature);
    if (!coreResistanceCoefficients.empty()) {
        circuitString += emit_core_ladder_spice(coreResistanceCoefficients, coil.get_functional_description().size());
    }

    // Mutual resistance: auxiliary winding network for cross-coupling losses (Hesterman 2020)
    // Only for multi-winding transformers when setting is enabled
    size_t numWindings = coil.get_functional_description().size();
    bool includeMutualResistance = settings.get_circuit_simulator_include_mutual_resistance();
    if (includeMutualResistance && numWindings >= 2) {
        auto mutualResistanceCoeffs = CircuitSimulatorExporter::calculate_mutual_resistance_coefficients(magnetic, temperature);
        if (!mutualResistanceCoeffs.empty()) {
            circuitString += emit_mutual_resistance_network_spice(mutualResistanceCoeffs, magnetizingInductance, numWindings);
        }
    }

    return headerString + "\n" + circuitString + "\n" + parametersString + "\n" + footerString;
}

std::string CircuitSimulatorExporterLtspiceModel::export_magnetic_as_symbol(Magnetic magnetic, std::optional<std::string> filePathOrFile) {
    std::string symbolString = "Version 4\n";
    symbolString += "SymbolType BLOCK\n";

    auto coil = magnetic.get_coil();

    int rectangleSemiWidth = 72;

    int leftSideSize = 16;
    int rightSideSize = 16;
    for (size_t index = 0; index < coil.get_functional_description().size(); index++) {
        if (coil.get_functional_description()[index].get_isolation_side() == IsolationSide::PRIMARY) {
            leftSideSize += 64;
        }
        else {
            rightSideSize += 64;
        }
    }

    int rectangleHeight = std::max(leftSideSize, rightSideSize);


    symbolString += "TEXT " + std::to_string(-rectangleSemiWidth + 8) + " " + std::to_string(-rectangleHeight / 2 + 8) + " Left 0 " + magnetic.get_reference() + "\n";
    symbolString += "TEXT " + std::to_string(-rectangleSemiWidth + 8) + " " + std::to_string(rectangleHeight / 2 - 8) + " Left 0 Made with OpenMagnetics\n";


    symbolString += "RECTANGLE Normal " + std::to_string(-rectangleSemiWidth) + " -" + std::to_string(rectangleHeight / 2) + " " + std::to_string(rectangleSemiWidth) + " " + std::to_string(rectangleHeight / 2) + "\n";
    symbolString += "SYMATTR Prefix X\n";
    symbolString += "SYMATTR Value " + fix_filename(magnetic.get_reference()) + "\n";
    symbolString += "SYMATTR ModelFile " + fix_filename(magnetic.get_reference()) + ".cir\n";

    int currentSpiceOrder = 1;
    int currentRectangleLeftSideHeight = -leftSideSize / 2 + 24;
    int currentRectangleRightSideHeight = -rightSideSize / 2 + 24;
    for (size_t index = 0; index < coil.get_functional_description().size(); index++) {

        if (coil.get_functional_description()[index].get_isolation_side() == IsolationSide::PRIMARY) {
            symbolString += "PIN " + std::to_string(-rectangleSemiWidth) + " " + std::to_string(currentRectangleLeftSideHeight) + " LEFT 8\n";
            currentRectangleLeftSideHeight += 32;
        }
        else{
            symbolString += "PIN " + std::to_string(rectangleSemiWidth) + " " + std::to_string(currentRectangleRightSideHeight) + " RIGHT 8\n";
            currentRectangleRightSideHeight += 32;
        }
        symbolString += "PINATTR PinName P" + std::to_string(index + 1) + "+\n";
        symbolString += "PINATTR SpiceOrder " + std::to_string(currentSpiceOrder) + "\n";
        currentSpiceOrder++;

        if (coil.get_functional_description()[index].get_isolation_side() == IsolationSide::PRIMARY) {
            symbolString += "PIN " + std::to_string(-rectangleSemiWidth) + " " + std::to_string(currentRectangleLeftSideHeight) + " LEFT 8\n";
            currentRectangleLeftSideHeight += 32;
        }
        else {
            symbolString += "PIN " + std::to_string(rectangleSemiWidth) + " " + std::to_string(currentRectangleRightSideHeight) + " RIGHT 8\n";
            currentRectangleRightSideHeight += 32;
        }
        symbolString += "PINATTR PinName P" + std::to_string(index + 1) + "-\n";
        symbolString += "PINATTR SpiceOrder " + std::to_string(currentSpiceOrder) + "\n";
        currentSpiceOrder++;
    }


    return symbolString;
}

void CircuitSimulationReader::process_line(std::string line, char separator) {
    std::stringstream ss(line);
    std::string token;

    if (_columns.size() == 0) {
        // Getting column names
        while(getline(ss, token, separator)) {
            if (int(token[0]) < 32) {
                continue;
            }
            CircuitSimulationSignal circuitSimulationSignal;
            token.erase(std::remove(token.begin(), token.end(), (char)13), token.end());
            token.erase(std::remove(token.begin(), token.end(), '\"'), token.end());
            circuitSimulationSignal.name = token;
            _columns.push_back(circuitSimulationSignal);
        }
    }
    else {
        size_t currentColumnIndex = 0;
        while(getline(ss, token, separator)) {
            if (int(token[0]) < 32) {
                continue;
            }
            _columns[currentColumnIndex].data.push_back(stod(token));
            currentColumnIndex++;
        }
    }
}

CircuitSimulationReader::CircuitSimulationReader(std::string filePathOrFile, bool forceFile) {
    char separator = '\0';
    std::string line;

    std::filesystem::path path = filePathOrFile;

    if (path.has_parent_path() && !forceFile) {
        if (!std::filesystem::exists(filePathOrFile)) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "File not found");
        }
        std::ifstream is(filePathOrFile);
        if (is.is_open()) {
            while(getline(is, line)) {
                if (separator == '\0') {
                    separator = guess_separator(line);
                }

                process_line(line, separator);
            }
            is.close();
        }
        else {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "File not found");
        }
    }
    else {
        std::stringstream ss(filePathOrFile);
        while(std::getline(ss, line, '\n')){
            if (separator == '\0') {
                separator = guess_separator(line);
            }

            process_line(line, separator);
        }
    }

    _time = find_time(_columns);
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

    for (size_t index = 1; index < data.size(); index++) {
        if (data[index - 1] >= data[index]) {
            return false;
        }
    }

    return true;
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
    std::vector<char> possibleSeparators = {',', ';', '\t'};

    for (auto separator : possibleSeparators) {
        std::stringstream ss(line);
        std::string token;
        size_t numberColumns = 0;
        while(getline(ss, token, separator)) {
            numberColumns++;
        }

        if (numberColumns >= 2) {
            return separator;
        }
    }

    throw InvalidInputException(ErrorCode::INVALID_INPUT, "No column separator found");
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
    auto& settings = Settings::GetInstance();
    Waveform waveform;
    waveform.set_data(signal.data);
    waveform.set_time(_time.data);

    auto waveformOnePeriod = get_one_period(waveform, frequency, sample);
    Waveform reconstructedWaveform = waveformOnePeriod;

    if (false) {
        double originalThreshold = settings.get_harmonic_amplitude_threshold();
        while (reconstructedWaveform.get_data().size() > 8192) {
            settings.set_harmonic_amplitude_threshold(settings.get_harmonic_amplitude_threshold() * 2);
            auto harmonics = Inputs::calculate_harmonics_data(waveformOnePeriod, frequency);
            settings.set_inputs_number_points_sampled_waveforms(2 * OpenMagnetics::round_up_size_to_power_of_2(harmonics.get_frequencies().back() / frequency));
            reconstructedWaveform = Inputs::reconstruct_signal(harmonics, frequency);
        }

        settings.set_harmonic_amplitude_threshold(originalThreshold);
        return reconstructedWaveform;
    }
    else {
        return waveformOnePeriod;
    }
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
    std::map<std::string, size_t> windingLabels = {
        {"pri", 0},
        {"sec", 1},
        {"aux", 2},
        {"ter", 2},
        {"a", 0},
        {"b", 1},
        {"c", 2},
        {"HV", 0},
        {"LV", 1},
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
                excitation.set_current(current);
            }
            if (column.windingIndex == windingIndex && column.type == DataType::MAGNETIZING_CURRENT) {
                auto waveform = extract_waveform(column, frequency);
                SignalDescriptor current;
                current.set_waveform(waveform);
                // Calculate processed data for magnetizing current
                auto currentProcessed = Inputs::calculate_processed_data(waveform, frequency, true);
                current.set_processed(currentProcessed);
                excitation.set_magnetizing_current(current);
            }
            if (column.windingIndex == windingIndex && column.type == DataType::VOLTAGE) {
                auto waveform = extract_waveform(column, frequency);
                SignalDescriptor voltage;
                voltage.set_waveform(waveform);
                // Calculate processed data for voltage
                auto voltageProcessed = Inputs::calculate_processed_data(waveform, frequency, true);
                voltage.set_processed(voltageProcessed);
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

    // Extract Steinmetz alpha from core material
    double steinmetzAlpha = 1.5;  // Default for typical ferrite
    try {
        auto core = magnetic.get_core();
        auto material = core.resolve_material();
        double refFreq = std::sqrt(opts.f0 * opts.f1);
        auto steinmetzDatum = CoreLossesModel::get_steinmetz_coefficients(material, refFreq);
        steinmetzAlpha = steinmetzDatum.get_alpha();
    } catch (...) {
        // Use default if material data unavailable
    }

    // alpha for fracpole: core resistance R_core ∝ f^(steinmetz_alpha - 1)
    opts.alpha = std::clamp(steinmetzAlpha - 1.0, 0.05, 0.95);

    // Choose profile: DD for alpha <= 0.5, FF for alpha > 0.5 (per Kundert recommendation)
    opts.profile = (opts.alpha <= 0.5) ? FracpoleProfile::DD : FracpoleProfile::FF;

    Curve2D coreResData = Sweeper().sweep_core_resistance_over_frequency(
        magnetic, opts.f0, opts.f1, 5, temperature);
    auto coreResVec = coreResData.get_y_points();
    auto freqVec    = coreResData.get_x_points();

    if (coreResVec.empty()) {
        // Return default network if no data
        return FractionalPole::generate(opts);
    }

    size_t midIdx = coreResVec.size() / 2;
    double R_ref = coreResVec[midIdx];
    double f_ref_actual = freqVec[midIdx];

    opts.coef = FractionalPole::fit_coef_from_reference(
        f_ref_actual, R_ref, opts.alpha, opts.f0, opts.f1,
        opts.lumpsPerDecade, opts.profile);

    return FractionalPole::generate(opts);
}

// ============================================================================
// NL5 circuit simulator exporter
// ============================================================================

static std::string make_nl5_key(const std::string& seed) {
    // Deterministic pseudo-UUID from seed hash (not cryptographic, just unique)
    std::hash<std::string> hasher;
    size_t h1 = hasher(seed);
    size_t h2 = hasher(seed + "_b");
    size_t h3 = hasher(seed + "_c");
    size_t h4 = hasher(seed + "_d");
    char buf[33];
    snprintf(buf, sizeof(buf),
             "%08X%08X%08X%08X",
             static_cast<uint32_t>(h1),
             static_cast<uint32_t>(h2),
             static_cast<uint32_t>(h3),
             static_cast<uint32_t>(h4));
    return std::string(buf);
}

// ---- NL5 Component builder ----
struct Nl5Cmp {
    std::string type;      // e.g. "R_R", "L_L", "C_C", "K_K", "G_G", "W_T2"
    int id;
    std::string name;
    int node0, node1;
    int node2 = -1, node3 = -1; // for 4-terminal (K, G)
    // Value parameter: for R=r, L=l, C=c, K=k
    std::string valParam;  // parameter name ("r","l","c","k")
    std::string valStr;    // value as string
    // Extra parameters
    std::map<std::string, std::string> extra;
    // Position
    int px = 0, py = 0, angle = 0;
    bool flip = false;
};

// ---- Engineering notation formatter for NL5 _txt attributes ----
static std::string nl5_to_eng(double v, int sigfigs = 4) {
    if (v == 0.0) return "0";
    bool neg = v < 0;
    double av = std::abs(v);
    std::string suffix; double scaled;
    if      (av >= 1e12) { scaled = av/1e12; suffix = "T"; }
    else if (av >= 1e9)  { scaled = av/1e9;  suffix = "G"; }
    else if (av >= 1e6)  { scaled = av/1e6;  suffix = "Meg"; }
    else if (av >= 1e3)  { scaled = av/1e3;  suffix = "k"; }
    else if (av >= 1.0)  { scaled = av;       suffix = ""; }
    else if (av >= 1e-3) { scaled = av*1e3;   suffix = "m"; }
    else if (av >= 1e-6) { scaled = av*1e6;   suffix = "u"; }
    else if (av >= 1e-9) { scaled = av*1e9;   suffix = "n"; }
    else if (av >= 1e-12){ scaled = av*1e12;  suffix = "p"; }
    else                 { scaled = av*1e15;  suffix = "f"; }
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*g", sigfigs, scaled);
    return (neg ? "-" : "") + std::string(buf) + suffix;
}

static std::string nl5_cmp_to_xml(const Nl5Cmp& cmp) {
    std::ostringstream ss;
    ss << "        <Cmp type=\"" << cmp.type << "\" id=\"" << cmp.id
       << "\" name=\"" << cmp.name << "\" descr=\"\" view=\"0\""
       << " node0=\"" << cmp.node0 << "\" node1=\"" << cmp.node1 << "\"";
    if (cmp.node2 >= 0) ss << " node2=\"" << cmp.node2 << "\"";
    if (cmp.node3 >= 0) ss << " node3=\"" << cmp.node3 << "\"";
    if (!cmp.valParam.empty()) {
        std::string mdl(1, static_cast<char>(std::toupper(static_cast<unsigned char>(cmp.valParam[0]))));
        ss << " model=\"" << mdl << "\"";
    }
    if (!cmp.valParam.empty()) {
        ss << " " << cmp.valParam << "=\"" << cmp.valStr << "\"";
        double dval = 0.0;
        try { dval = std::stod(cmp.valStr); } catch (...) {}
        ss << " " << cmp.valParam << "_txt=\"" << nl5_to_eng(dval) << "\"";
    }
    for (auto& [k, v] : cmp.extra) ss << " " << k << "=\"" << v << "\"";
    ss << ">\n";
    ss << "          <Pos x=\"" << cmp.px << "\" y=\"" << cmp.py
       << "\" angle=\"" << cmp.angle << "\" flip=\"" << (cmp.flip ? "1" : "0") << "\" />\n";
    ss << "        </Cmp>\n";
    return ss.str();
}

static std::string nl5_wire_to_xml(int id, int node0, int node1, int x0, int y0, int x1, int y1, const std::string& name="") {
    std::ostringstream ss;
    ss << "        <Cmp type=\"W_T2\" id=\"" << id << "\" name=\"" << (name.empty() ? "W" + std::to_string(id) : name)
       << "\" descr=\"\" view=\"0\" node0=\"" << node0 << "\" node1=\"" << node1 << "\">\n";
    ss << "          <Pos x=\"" << x0 << "\" y=\"" << y0 << "\" angle=\"0\" flip=\"0\" />\n";
    ss << "          <Point x=\"" << x1 << "\" y=\"" << y1 << "\" />\n";
    ss << "        </Cmp>\n";
    return ss.str();
}

// ============================================================================
// CircuitSimulatorExporterNl5Model::export_magnetic_as_subcircuit
// ============================================================================
// Uses NL5 Winding (W) components as documented in NL5 Reference:
// - Each winding is a W_Winding component with turns ratio
// - All winding "core" pins connect to a common core node
// - Magnetizing inductance connects from core node to ground
// - Leakage is modeled via series inductors
// - DC resistance and AC losses via series R/L/C networks
// ============================================================================
std::string CircuitSimulatorExporterNl5Model::export_magnetic_as_subcircuit(
    Magnetic magnetic,
    double frequency,
    double temperature,
    std::optional<std::string> filePathOrFile,
    CircuitSimulatorExporterCurveFittingModes mode)
{
    auto coil = magnetic.get_coil();
    size_t numWindings = coil.get_functional_description().size();

    // Resolve mode
    auto resolvedMode = CircuitSimulatorExporter::resolve_curve_fitting_mode(mode);

    // Precompute ladder coefficients (if needed)
    std::vector<std::vector<double>> acCoeffs;
    std::vector<FractionalPoleNetwork> fracNets;
    std::optional<FractionalPoleNetwork> coreFracNet;

    if (resolvedMode == CircuitSimulatorExporterCurveFittingModes::FRACPOLE) {
        fracNets = CircuitSimulatorExporter::calculate_fracpole_networks_per_winding(magnetic, temperature);
        try { coreFracNet = CircuitSimulatorExporter::calculate_core_fracpole_network(magnetic, temperature); }
        catch (...) {}
    } else if (resolvedMode != CircuitSimulatorExporterCurveFittingModes::ANALYTICAL) {
        acCoeffs = CircuitSimulatorExporter::calculate_ac_resistance_coefficients_per_winding(
            magnetic, temperature, CircuitSimulatorExporterCurveFittingModes::LADDER);
    }

    auto coreCoeffs = CircuitSimulatorExporter::calculate_core_resistance_coefficients(magnetic, temperature);

    double magnetizingInductance = resolve_dimensional_values(
        MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(magnetic)
            .get_magnetizing_inductance());
    auto leakageInductances = LeakageInductance()
        .calculate_leakage_inductance(magnetic, Defaults().measurementFrequency)
        .get_leakage_inductance_per_winding();

    // Get primary turns for reference
    double primaryTurns = static_cast<double>(coil.get_functional_description()[0].get_number_turns());
    
    // Check if saturation modeling is enabled
    auto& settings = Settings::GetInstance();
    bool includeSaturation = settings.get_circuit_simulator_include_saturation();
    SaturationParameters satParams;
    if (includeSaturation) {
        satParams = get_saturation_parameters(magnetic, temperature);
        // Saturation only supported for single winding (inductor) - coupled windings can't use saturating inductors
        if (satParams.valid && numWindings > 1) {
            // For transformers, fall back to linear model with warning
            // (NL5 XML comment not easily added, so just disable saturation)
            includeSaturation = false;
        }
    }

    // ---- NL5 XML generation ----
    // Node numbering scheme:
    //   0          = ground
    //   1..N       = winding positive terminal Pi+  (external)
    //   N+1..2N    = winding negative terminal Pi-  (external)
    //   50         = core node (common for all windings via ideal transformer)
    //   100+i      = internal node after Rdc (before ladder/fracpole), winding i
    //   200+i      = internal node after ladder/fracpole (before leakage L), winding i
    //   250+i      = internal node after leakage (connects to winding), winding i
    //   300+i*20+k = internal nodes within fracpole/ladder stages for winding i, stage k
    //   900+k      = internal nodes for core loss network

    int nextId = 1;   // component id counter
    std::string cmps; // accumulated <Cmp> XML

    // Core node where all windings connect via ideal transformers
    const int node_core = 50;

    auto add_R = [&](const std::string& name, int n0, int n1, double val, int px, int py, int angle=0) {
        Nl5Cmp c; c.type="R_R"; c.id=nextId++; c.name=name;
        c.node0=n0; c.node1=n1; c.valParam="r"; c.valStr=to_string(val,12);
        c.px=px; c.py=py; c.angle=angle;
        cmps += nl5_cmp_to_xml(c);
    };
    auto add_L = [&](const std::string& name, int n0, int n1, double val, int px, int py, int angle=0) {
        Nl5Cmp c; c.type="L_L"; c.id=nextId++; c.name=name;
        c.node0=n0; c.node1=n1; c.valParam="l"; c.valStr=to_string(val,12);
        c.px=px; c.py=py; c.angle=angle;
        cmps += nl5_cmp_to_xml(c);
    };
    auto add_C = [&](const std::string& name, int n0, int n1, double val, int px, int py, int angle=0) {
        Nl5Cmp c; c.type="C_C"; c.id=nextId++; c.name=name;
        c.node0=n0; c.node1=n1; c.valParam="c"; c.valStr=to_string(val,12);
        c.px=px; c.py=py; c.angle=angle;
        cmps += nl5_cmp_to_xml(c);
    };
    // NL5 saturating inductor using L_Lsat model
    // Parameters: l=initial inductance, Lsat=saturated inductance, Iknee=knee current
    auto add_Lsat = [&](const std::string& name, int n0, int n1, double L0, double Lsat, double Iknee,
                        int px, int py, int angle=0) {
        Nl5Cmp c; c.type="L_Lsat"; c.id=nextId++; c.name=name;
        c.node0=n0; c.node1=n1; c.valParam="l"; c.valStr=to_string(L0, 12);
        // Lsat is the saturated (minimum) inductance
        c.extra["Lsat"] = to_string(Lsat, 12);
        c.extra["Lsat_txt"] = nl5_to_eng(Lsat);
        // Iknee is the current at knee point (where saturation starts)
        c.extra["Iknee"] = to_string(Iknee, 6);
        c.extra["Iknee_txt"] = nl5_to_eng(Iknee);
        c.extra["model"] = "Lsat";
        c.px=px; c.py=py; c.angle=angle;
        cmps += nl5_cmp_to_xml(c);
    };
    // NL5 Winding component (W_Winding): ideal transformer with turns ratio
    // node0 = winding +, node1 = winding -, node2 = core
    auto add_Winding = [&](const std::string& name, int n_plus, int n_minus, int n_core, 
                           double turns, int px, int py) {
        Nl5Cmp c; c.type="W_Winding"; c.id=nextId++; c.name=name;
        c.node0=n_plus; c.node1=n_minus; c.node2=n_core;
        c.valParam="n"; c.valStr=to_string(turns, 6);
        c.extra["model"] = "Winding";
        c.px=px; c.py=py;
        cmps += nl5_cmp_to_xml(c);
    };

    // Layout: each winding is 80 units wide, stacked vertically
    for (size_t wi = 0; wi < numWindings; ++wi) {
        int base_y = static_cast<int>(wi) * 80;
        int is = static_cast<int>(wi) + 1;
        std::string ws = std::to_string(is);

        int node_Pplus  = static_cast<int>(wi) + 1;           // Pi+  external
        int node_Pminus = static_cast<int>(numWindings) + is;  // Pi-  external
        int node_after_rdc  = 100 + is;                        // after Rdc
        int node_after_net  = 200 + is;                        // after ladder/fracpole
        int node_after_leak = 250 + is;                        // after leakage inductance

        double Rdc = WindingLosses::calculate_effective_resistance_of_winding(magnetic, wi, 0.1, temperature);
        double numTurns = static_cast<double>(coil.get_functional_description()[wi].get_number_turns());
        
        // Get leakage inductance for this winding (referred to primary)
        double leakageL = 0.0;
        if (wi > 0 && wi <= leakageInductances.size()) {
            leakageL = resolve_dimensional_values(leakageInductances[wi-1]);
            // Scale leakage to this winding's turns
            double turnsRatio = numTurns / primaryTurns;
            leakageL *= turnsRatio * turnsRatio;
        }

        // Add port wire labels (external terminals)
        cmps += nl5_wire_to_xml(nextId++, node_Pplus, node_Pplus, 0, base_y, -20, base_y, "P"+ws+"+");
        cmps += nl5_wire_to_xml(nextId++, node_Pminus, node_Pminus, 180, base_y+40, 200, base_y+40, "P"+ws+"-");

        // 1. Rdc (DC resistance)
        add_R("Rdc"+ws, node_Pplus, node_after_rdc, Rdc, 20, base_y);

        int node_last = node_after_rdc;

        // 2. AC loss network (FRACPOLE or LADDER)
        if (resolvedMode == CircuitSimulatorExporterCurveFittingModes::FRACPOLE && wi < fracNets.size()) {
            // FRACPOLE: R/C stages
            const auto& net = fracNets[wi];
            for (size_t k = 0; k < net.stages.size(); ++k) {
                int node_stage = 300 + is * 20 + static_cast<int>(k);
                int px_R = 40 + static_cast<int>(k) * 16;
                // Series R from last node to stage node
                add_R("Rfp"+ws+"_"+std::to_string(k), node_last, node_stage,
                      net.stages[k].R * net.opts.coef, px_R, base_y);
                // Shunt C from stage node to ground
                add_C("Cfp"+ws+"_"+std::to_string(k), node_stage, 0,
                      net.stages[k].C / net.opts.coef, px_R+4, base_y+20, 90);
                node_last = node_stage;
            }
            if (net.Cinf > 0.0) {
                add_C("Cfpinf"+ws, node_last, 0, net.Cinf / net.opts.coef,
                      40 + static_cast<int>(net.stages.size())*16, base_y+20, 90);
            }
            // Wire to next stage
            cmps += nl5_wire_to_xml(nextId++, node_last, node_after_net,
                  40 + static_cast<int>(net.stages.size())*16 + 4, base_y, 
                  60 + static_cast<int>(net.stages.size())*16, base_y);
            node_last = node_after_net;
        } else if (resolvedMode == CircuitSimulatorExporterCurveFittingModes::LADDER
                   && !acCoeffs.empty() && wi < acCoeffs.size()) {
            // LADDER: R/L stages
            const auto& c = acCoeffs[wi];
            size_t stages = c.size() / 2;
            for (size_t k = 0; k < stages; ++k) {
                int node_stage = 300 + is * 20 + static_cast<int>(k);
                int px = 40 + static_cast<int>(k) * 16;
                add_L("Ll"+ws+"_"+std::to_string(k), node_last, node_stage, c[2*k], px, base_y);
                add_R("Rl"+ws+"_"+std::to_string(k), node_stage, 0, c[2*k+1], px+4, base_y+20, 90);
                node_last = node_stage;
            }
            // Wire to next stage
            cmps += nl5_wire_to_xml(nextId++, node_last, node_after_net,
                  40 + static_cast<int>(stages)*16, base_y, 60 + static_cast<int>(stages)*16, base_y);
            node_last = node_after_net;
        } else {
            // ANALYTICAL or fallback: just wire node_after_rdc to node_after_net
            cmps += nl5_wire_to_xml(nextId++, node_after_rdc, node_after_net,
                                    40, base_y, 80, base_y);
            node_last = node_after_net;
        }

        // 3. Leakage inductance (series with winding)
        if (leakageL > 0.0 && wi > 0) {
            add_L("Lleak"+ws, node_last, node_after_leak, leakageL, 100, base_y);
            node_last = node_after_leak;
        } else {
            // No leakage for primary or if leakage is zero
            node_after_leak = node_last;
        }

        // 4. Winding (ideal transformer to core)
        // NL5 Winding: node0=+, node1=-, node2=core
        add_Winding("W"+ws, node_after_leak, node_Pminus, node_core, numTurns, 130, base_y);
    }

    // 5. Magnetizing inductance (from core to ground)
    // This represents the magnetic energy storage in the core
    if (includeSaturation && satParams.valid) {
        // Use saturating inductor model
        const double mu0 = 4e-7 * M_PI;
        double Isat = satParams.Isat();
        // Calculate gap contribution for saturated inductance
        double effectiveMuR = satParams.Lmag / (mu0 * satParams.primaryTurns * satParams.primaryTurns * satParams.Ae / satParams.le);
        double gapFactor = std::min(0.99, 1.0 / effectiveMuR);
        double Lsat = satParams.Lmag * gapFactor;  // Saturated inductance (just the gap)
        // Ensure reasonable minimum values
        if (Isat < 1e-6) Isat = 1e-6;
        if (Lsat < 1e-12) Lsat = magnetizingInductance * 0.01;
        // NL5 Lsat model: L0 is initial inductance, Lsat is saturated inductance, Iknee is knee current
        add_Lsat("Lmag", node_core, 0, magnetizingInductance, Lsat, Isat,
                 160, static_cast<int>(numWindings) * 40);
    } else {
        add_L("Lmag", node_core, 0, magnetizingInductance, 160, static_cast<int>(numWindings) * 40);
    }

    // 6. Core loss network in parallel with Lmag using R-L ladder
    // The R-L ladder provides impedance that INCREASES with frequency, matching core loss behavior.
    // coreCoeffs already calculated earlier in this function
    if (!coreCoeffs.empty()) {
        // Build R-L ladder: series of (L || R) stages from node_core to ground
        // Z = L1 || (R1 + L2 || (R2 + L3 || R3))
        int node_prev = 0;  // ground
        int stageCount = static_cast<int>(coreCoeffs.size() / 2);
        for (int k = stageCount - 1; k >= 0; --k) {
            double R = coreCoeffs[k * 2];
            double L = coreCoeffs[k * 2 + 1];
            // Validate values
            if (R <= 0 || L <= 0 || R > 1e6 || L > 1) {
                continue;  // Skip invalid stage
            }
            int node_stage = 900 + k;
            int y_offset = static_cast<int>(numWindings) * 40 + k * 16;
            // R in series, L in parallel
            add_R("Rcore" + std::to_string(k), node_stage, node_prev, R, 180, y_offset);
            int node_next = (k == 0) ? node_core : (900 + k - 1);
            add_L("Lcore" + std::to_string(k), node_next, node_prev, L, 184, y_offset + 8);
            node_prev = node_stage;
        }
    }

    // ---- Assemble full XML ----
    std::string ref = fix_filename(magnetic.get_reference());
    std::string key = make_nl5_key(ref);

    std::ostringstream xml;
    xml << "<?xml version=\"1.0\"?>\n";
    xml << "<NL5>\n";
    xml << "  <c>Generated by OpenMagnetics - do not edit manually</c>\n";
    xml << "  <c>Component: " << magnetic.get_reference() << "</c>\n";
    xml << "  <Version ver=\"3\" rev=\"18\" core=\"78\" build=\"99\" />\n";
    xml << "  <Doc>\n";
    xml << "    <Cir mode=\"2\">\n";
    xml << "      <Cmps Key=\"" << key << "\">\n";
    xml << cmps;
    xml << "      </Cmps>\n";
    xml << "      <Scopes />\n";
    xml << "      <Scripts />\n";
    xml << "    </Cir>\n";
    xml << "  </Doc>\n";
    xml << "</NL5>\n";

    std::string result = xml.str();

    if (filePathOrFile.has_value()) {
        std::ofstream outFile(filePathOrFile.value());
        outFile << result;
    }
    return result;
}

} // namespace OpenMagnetics

