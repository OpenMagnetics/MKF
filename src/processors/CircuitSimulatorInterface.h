#pragma once
#include "json.hpp"
#include "MAS.hpp"
#include "constructive_models/Magnetic.h"
#include "support/Utils.h"
#include <random>
#include <chrono>
#include "support/Exceptions.h"

using ordered_json = nlohmann::ordered_json;

using namespace MAS;


// ============================================================================
// SECTION: NonlinearLeastSquares
// Eigen-based Levenberg-Marquardt solver replacing levmar (no Fortran/LAPACK).
// ============================================================================

#include <Eigen/Dense>
#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/NumericalDiff>
#include <vector>
#include <cmath>
#include <limits>

namespace OpenMagnetics {

struct LMFunctorWrapper : Eigen::DenseFunctor<double> {
    using LevmarFunc = void (*)(double *p, double *x, int m, int n, void *data);

    LevmarFunc _func;
    const double* _measurements;
    void* _data;

    LMFunctorWrapper(int nParams, int nResiduals, LevmarFunc func,
                     const double* measurements, void* data)
        : Eigen::DenseFunctor<double>(nParams, nResiduals)
        , _func(func)
        , _measurements(measurements)
        , _data(data) {}

    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
        std::vector<double> modelOutput(static_cast<size_t>(values()));
        std::vector<double> paramsCopy(x.data(), x.data() + inputs());

        _func(paramsCopy.data(), modelOutput.data(), static_cast<int>(inputs()), static_cast<int>(values()), _data);

        for (int i = 0; i < values(); ++i) {
            fvec[i] = modelOutput[static_cast<size_t>(i)] - _measurements[i];
        }
        return 0;
    }
};

inline int eigen_levmar_dif(
    void (*func)(double *p, double *x, int m, int n, void *data),
    double *p, double *x, int m, int n, int itmax,
    double *opts, double *info,
    [[maybe_unused]] double *work,
    [[maybe_unused]] double *covar,
    void *data)
{
    LMFunctorWrapper functor(m, n, func, x, data);
    Eigen::NumericalDiff<LMFunctorWrapper> numDiff(functor);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<LMFunctorWrapper>> lm(numDiff);

    // Set parameters using the public member directly
    lm.setMaxfev(itmax * (m + 1));
    if (opts != nullptr) {
        lm.setFtol(opts[1]);
        lm.setXtol(opts[2]);
        lm.setGtol(opts[3]);
    } else {
        lm.setFtol(1e-25);
        lm.setXtol(1e-25);
        lm.setGtol(1e-25);
    }

    Eigen::VectorXd params(m);
    for (int i = 0; i < m; ++i) {
        params[i] = p[i];
    }

    double initialNorm = 0.0;
    if (info != nullptr) {
        Eigen::VectorXd initialResiduals(n);
        functor(params, initialResiduals);
        initialNorm = initialResiduals.squaredNorm();
    }

    auto status = lm.minimize(params);

    for (int i = 0; i < m; ++i) {
        p[i] = params[i];
    }

    if (info != nullptr) {
        info[0] = initialNorm;
        info[1] = lm.fvec().squaredNorm();
        info[2] = 0;
        info[3] = 0;
        info[4] = 0;
        info[5] = static_cast<double>(lm.iterations());
        info[6] = static_cast<double>(static_cast<int>(status));
        info[7] = static_cast<double>(lm.nfev());
        info[8] = static_cast<double>(lm.njev());
        info[9] = 0;
    }

    return static_cast<int>(lm.iterations());
}

} // namespace OpenMagnetics


// ============================================================================
// SECTION: FractionalPole
// Analytical fractional-pole R/C generator for winding and core loss modeling.
// ============================================================================

#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <numbers>
#include <sstream>
#include <iomanip>

namespace OpenMagnetics {

enum class FracpoleProfile {
    DD, DF, FD, FF
};

inline FracpoleProfile fracpole_profile_from_string(const std::string& s) {
    if (s == "dd" || s == "DD") return FracpoleProfile::DD;
    if (s == "df" || s == "DF") return FracpoleProfile::DF;
    if (s == "fd" || s == "FD") return FracpoleProfile::FD;
    if (s == "ff" || s == "FF") return FracpoleProfile::FF;
    throw std::invalid_argument("Unknown FracpoleProfile: " + s + " (use dd/df/fd/ff)");
}

inline std::string fracpole_profile_to_string(FracpoleProfile p) {
    switch (p) {
        case FracpoleProfile::DD: return "dd";
        case FracpoleProfile::DF: return "df";
        case FracpoleProfile::FD: return "fd";
        case FracpoleProfile::FF: return "ff";
        default: return "dd";
    }
}

struct FractionalPoleOptions {
    double f0             = 0.1;
    double f1             = 1e7;
    double alpha          = 0.5;
    double lumpsPerDecade = 1.5;
    FracpoleProfile profile = FracpoleProfile::DD;
    double coef           = 1.0;
};

struct FractionalPoleStage {
    double R;
    double C;
};

struct FractionalPoleNetwork {
    std::vector<FractionalPoleStage> stages;
    double Cinf = 0.0;
    double Rinf = 0.0;
    FractionalPoleOptions opts;
};

class FractionalPole {
public:
    static FractionalPoleNetwork generate(const FractionalPoleOptions& opts) {
        const double f0     = opts.f0;
        const double f1     = opts.f1;
        const double alpha  = opts.alpha;
        const double lpd    = opts.lumpsPerDecade;
        const FracpoleProfile profile = opts.profile;
        // coef is stored in opts but NOT used during generation.
        // R/C values are base (unscaled) values per Kundert convention.
        // coef is applied in evaluate_*() and to_spice_*() functions.

        if (f0 <= 0 || f1 <= f0)
            throw std::invalid_argument("FractionalPole: require 0 < f0 < f1");
        if (alpha <= 0.0 || alpha >= 1.0)
            throw std::invalid_argument("FractionalPole: alpha must be in (0,1)");
        if (lpd <= 0.0)
            throw std::invalid_argument("FractionalPole: lumpsPerDecade must be > 0");

        const double decades = std::log10(f1 / f0);
        const int N = static_cast<int>(std::ceil(lpd * decades));
        const double r = std::pow(f1 / f0, 1.0 / N);
        const double r_pz = std::pow(r, alpha);

        double f_zero_start = 0.0;
        switch (profile) {
            case FracpoleProfile::DD:
            case FracpoleProfile::DF:
                f_zero_start = f0;
                break;
            case FracpoleProfile::FD:
            case FracpoleProfile::FF:
                f_zero_start = f0 * r_pz;
                break;
        }

        FractionalPoleNetwork net;
        net.opts = opts;

        const double twoPi = 2.0 * std::numbers::pi;

        // Compute BASE R/C values (without coef scaling).
        // Kundert convention: actual R = base_R * coef, actual C = base_C / coef
        // This ensures coef acts as a pure impedance scaling factor.
        for (int i = 0; i < N; ++i) {
            double fz = f_zero_start * std::pow(r, i);
            double fp_stage = fz * r_pz;
            double fc = std::sqrt(fz * fp_stage);

            double stage_R = std::sqrt(r_pz) / (twoPi * fc);
            double stage_C = 1.0 / (std::sqrt(r_pz) * twoPi * fc);

            net.stages.push_back({stage_R, stage_C});
        }

        // Terminal elements (also base values, without coef)
        switch (profile) {
            case FracpoleProfile::DD: {
                double last_fz = f_zero_start * std::pow(r, N - 1);
                double last_fp = last_fz * r_pz;
                if (!net.stages.empty()) {
                    double Rlast = net.stages.back().R;
                    net.Cinf = 1.0 / (twoPi * last_fp * Rlast);
                }
                break;
            }
            case FracpoleProfile::FF: {
                double first_fp = f0;
                if (!net.stages.empty()) {
                    double Cfirst = net.stages.front().C;
                    net.Rinf = 1.0 / (twoPi * first_fp * Cfirst);
                }
                break;
            }
            default:
                break;
        }

        return net;
    }

    static FractionalPoleNetwork generate(double f0, double f1, double alpha,
                                          double lumpsPerDecade = 1.5,
                                          FracpoleProfile profile = FracpoleProfile::DD,
                                          double coef = 1.0) {
        FractionalPoleOptions opts;
        opts.f0 = f0; opts.f1 = f1; opts.alpha = alpha;
        opts.lumpsPerDecade = lumpsPerDecade;
        opts.profile = profile; opts.coef = coef;
        return generate(opts);
    }

    // -------------------------------------------------------------------------
    // evaluate_impedance()
    //
    // Evaluates |Z| of the lumped network at frequency f.
    // Kundert convention: R_actual = R_base * coef, C_actual = C_base / coef
    // so Y_R = 1/(R_base * coef), Y_C = jw * C_base / coef
    // -------------------------------------------------------------------------
    static double evaluate_impedance(const FractionalPoleNetwork& net, double f) {
        const double w = 2.0 * std::numbers::pi * f;
        const double coef = net.opts.coef;
        std::complex<double> Y_total(0.0, 0.0);
        for (const auto& s : net.stages) {
            // Each stage: R in series with C, both in parallel with the port
            // Admittance of series RC: Y = 1 / (R*coef + 1/(jw*C/coef))
            //   = jw*C/coef / (1 + jw*R*coef*C/coef)
            //   = jw*C/coef / (1 + jw*R*C)
            double R_act = s.R * coef;
            double C_act = s.C / coef;
            std::complex<double> Z_branch(R_act, -1.0 / (w * C_act));
            Y_total += 1.0 / Z_branch;
        }
        if (net.Cinf > 0.0)
            Y_total += std::complex<double>(0.0, w * net.Cinf / coef);
        if (net.Rinf > 0.0)
            Y_total += std::complex<double>(1.0 / (net.Rinf * coef), 0.0);

        return std::abs(1.0 / Y_total);
    }

    // -------------------------------------------------------------------------
    // evaluate_resistance()
    //
    // Returns Re(Z) = effective AC resistance at frequency f.
    // Kundert convention: R_actual = R_base * coef, C_actual = C_base / coef
    // -------------------------------------------------------------------------
    static double evaluate_resistance(const FractionalPoleNetwork& net, double f) {
        const double w = 2.0 * std::numbers::pi * f;
        const double coef = net.opts.coef;
        std::complex<double> Y_total(0.0, 0.0);
        for (const auto& s : net.stages) {
            double R_act = s.R * coef;
            double C_act = s.C / coef;
            std::complex<double> Z_branch(R_act, -1.0 / (w * C_act));
            Y_total += 1.0 / Z_branch;
        }
        if (net.Cinf > 0.0)
            Y_total += std::complex<double>(0.0, w * net.Cinf / coef);
        if (net.Rinf > 0.0)
            Y_total += std::complex<double>(1.0 / (net.Rinf * coef), 0.0);

        std::complex<double> Z = 1.0 / Y_total;
        return Z.real();
    }

    static double fit_coef_from_reference(double f_ref, double R_ref,
                                          double alpha,
                                          double f0, double f1,
                                          double lumpsPerDecade = 1.5,
                                          FracpoleProfile profile = FracpoleProfile::DD) {
        FractionalPoleNetwork unit_net = generate(f0, f1, alpha, lumpsPerDecade, profile, 1.0);
        double R_unit = evaluate_resistance(unit_net, f_ref);
        if (R_unit < 1e-30) return 1.0;
        return R_ref / R_unit;
    }

    static std::string to_spice_subcircuit(const FractionalPoleNetwork& net,
                                           const std::string& name,
                                           const std::string& portP = "p",
                                           const std::string& portN = "n") {
        std::ostringstream ss;
        ss << std::scientific << std::setprecision(6);
        ss << "* Fractional Impedance Pole (fracpole)\n";
        ss << "* Z ~ s^(-" << net.opts.alpha << "), valid from "
           << net.opts.f0 << " Hz to " << net.opts.f1 << " Hz\n";
        ss << "* Profile: " << fracpole_profile_to_string(net.opts.profile)
           << ", " << net.stages.size() << " stages\n";
        ss << ".subckt " << name << " " << portP << " " << portN << "\n";
        ss << ".param coef=" << net.opts.coef << "\n";

        for (size_t i = 0; i < net.stages.size(); ++i) {
            const auto& s = net.stages[i];
            std::string node = std::to_string(i + 1);
            ss << "R" << (i + 1) << " " << portP << " " << node
               << " r=" << s.R << "*coef\n";
            ss << "C" << (i + 1) << " " << node << " " << portN
               << " c=" << s.C << "/coef\n";
        }

        if (net.Cinf > 0.0) {
            ss << "Cinf " << portP << " " << portN
               << " c=" << net.Cinf << "/coef\n";
        }
        if (net.Rinf > 0.0) {
            ss << "Rinf " << portP << " " << portN
               << " r=" << net.Rinf << "*coef\n";
        }

        ss << ".ends " << name << "\n";
        return ss.str();
    }

    static std::string to_spice_skin_effect_subcircuit(
        const FractionalPoleNetwork& fracNet,
        const std::string& skinName,
        const std::string& fracpoleName,
        double s_param,
        double scaling = 1e6,
        double gmin = 1e-12)
    {
        std::ostringstream ss;
        ss << std::scientific << std::setprecision(6);
        ss << "* Skin effect subcircuit using fracpole + gyrator\n";
        ss << ".subckt " << skinName << " 1 2\n";
        ss << ".param s=" << s_param << "\n";
        ss << ".param scaling=" << scaling << "\n";
        ss << ".param gmin=" << gmin << "\n";
        ss << "* Gyrator: converts impedance fracpole -> admittance fracpole\n";
        ss << "G1 1 2 3 0 " << "cur='v(3,0)/sqrt(scaling)'\n";
        ss << "G2 0 3 2 1 " << "cur='v(2,1)/sqrt(scaling)'\n";
        ss << "Rg1 1 2 r='1/gmin'\n";
        ss << "Rg2 3 0 r='1/gmin'\n";
        ss << "X_fp 3 0 " << fracpoleName
           << " coef='scaling/s'\n";
        ss << ".ends " << skinName << "\n\n";

        ss << to_spice_subcircuit(fracNet, fracpoleName, "p", "n");
        return ss.str();
    }

    static int stage_count_for_accuracy(double f0, double f1,
                                        double lumpsPerDecade = 1.5) {
        return static_cast<int>(std::ceil(lumpsPerDecade * std::log10(f1 / f0)));
    }
};

} // namespace OpenMagnetics

namespace OpenMagnetics {


enum class CircuitSimulatorExporterModels : int {
    SIMBA,
    NGSPICE,
    LTSPICE,
    NL5,
    PLECS
};

void from_json(const json & j, CircuitSimulatorExporterModels & x);
void to_json(json & j, const CircuitSimulatorExporterModels & x);

inline void from_json(const json & j, CircuitSimulatorExporterModels & x) {
    if (j == "SIMBA") x = CircuitSimulatorExporterModels::SIMBA;
    else if (j == "NgSpice") x = CircuitSimulatorExporterModels::NGSPICE;
    else if (j == "LtSpice") x = CircuitSimulatorExporterModels::LTSPICE;
    else if (j == "NL5") x = CircuitSimulatorExporterModels::NL5;
    else if (j == "PLECS") x = CircuitSimulatorExporterModels::PLECS;
    else { throw std::runtime_error("Input JSON does not conform to schema!"); }
}

inline void to_json(json & j, const CircuitSimulatorExporterModels & x) {
    switch (x) {
        case CircuitSimulatorExporterModels::SIMBA: j = "SIMBA"; break;
        case CircuitSimulatorExporterModels::NGSPICE: j = "NgSpice"; break;
        case CircuitSimulatorExporterModels::LTSPICE: j = "LtSpice"; break;
        case CircuitSimulatorExporterModels::NL5: j = "NL5"; break;
        case CircuitSimulatorExporterModels::PLECS: j = "PLECS"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"[object Object]\": " + std::to_string(static_cast<int>(x)));
    }
}


enum class CircuitSimulatorExporterCurveFittingModes : int {
    ANALYTICAL,
    LADDER,
    FRACPOLE,
    AUTO
};

// Core loss ladder topology selection
// RIDLEY: RL, RL, RL nested stages (Ridley model)
// ROSANO: R || (R+L) || (R+L+C) parallel branches (Rosano model)
enum class CoreLossTopology : int {
    RIDLEY = 0,
    ROSANO = 1
};

class CircuitSimulatorExporterModel {
    private:
    protected:
        CircuitSimulatorExporterCurveFittingModes _acResistanceMode = CircuitSimulatorExporterCurveFittingModes::LADDER;
    public:
        std::string programName = "Default";
        virtual std::string export_magnetic_as_symbol(Magnetic magnetic, std::optional<std::string> filePathOrFile = std::nullopt) = 0;
        virtual std::string export_magnetic_as_subcircuit(Magnetic magnetic, double frequency = defaults.measurementFrequency, double temperature = defaults.ambientTemperature, std::optional<std::string> filePathOrFile = std::nullopt, CircuitSimulatorExporterCurveFittingModes mode=CircuitSimulatorExporterCurveFittingModes::LADDER) = 0;
        static std::shared_ptr<CircuitSimulatorExporterModel> factory(CircuitSimulatorExporterModels programName);

};

class CircuitSimulatorExporter {
    private:
    protected:
        std::shared_ptr<CircuitSimulatorExporterModel> _model;
    public:
        CircuitSimulatorExporter(CircuitSimulatorExporterModels program);
        CircuitSimulatorExporter();

        static double analytical_model(double x[], double frequency);
        static void analytical_func(double *p, double *x, int m, int n, void *data);
        static double ladder_model(double x[], double frequency, double dcResistance);
        static void ladder_func(double *p, double *x, int m, int n, void *data);
        static double core_ladder_model(double x[], double frequency, double dcResistance);
        static void core_ladder_func(double *p, double *x, int m, int n, void *data);
        static double core_rosano_model(double x[], double frequency);
        static void core_rosano_func(double *p, double *x, int m, int n, void *data);

        static double fracpole_skin_model(double x[], double frequency);
        static void fracpole_skin_func(double *p, double *x, int m, int n, void *data);
        static double fracpole_model(double x[], int numCoeffs, double frequency, double dcResistance);

        static std::vector<std::vector<double>> calculate_ac_resistance_coefficients_per_winding(Magnetic magnetic, double temperature, CircuitSimulatorExporterCurveFittingModes mode = CircuitSimulatorExporterCurveFittingModes::LADDER);
        
        /**
         * @brief Structure to hold mutual resistance ladder coefficients for a winding pair.
         * 
         * Based on Hesterman (2020) "Mutual Resistance" paper.
         * The mutual resistance R_ij represents the cross-coupling loss between windings.
         * 
         * Power dissipation: P = R11*I1² + R22*I2² + 2*R12*I1*I2*cos(θ)
         * where θ is the phase angle between currents.
         * 
         * For transformers, currents are typically ~180° out of phase, so:
         * - Positive R12: reduces total losses (cancellation)
         * - Negative R12: increases total losses (reinforcement)
         */
        struct MutualResistanceCoefficients {
            size_t windingIndex1;           // First winding index
            size_t windingIndex2;           // Second winding index
            std::vector<double> coefficients;  // R-L ladder coefficients [R1, L1, R2, L2, ...]
            double dcMutualResistance;      // DC value of mutual resistance
        };
        
        /**
         * @brief Calculate mutual resistance ladder coefficients for all winding pairs.
         * 
         * Returns a vector of MutualResistanceCoefficients, one for each unique pair.
         * For N windings, returns N*(N-1)/2 coefficient sets.
         * 
         * @param magnetic The magnetic component
         * @param temperature Operating temperature in Celsius
         * @return Vector of mutual resistance coefficients for each winding pair
         */
        static std::vector<MutualResistanceCoefficients> calculate_mutual_resistance_coefficients(
            Magnetic magnetic, double temperature);
        static std::vector<FractionalPoleNetwork> calculate_fracpole_networks_per_winding(Magnetic magnetic, double temperature, FractionalPoleOptions opts={});
        static FractionalPoleNetwork calculate_core_fracpole_network(Magnetic magnetic, double temperature, FractionalPoleOptions opts={});
        static CircuitSimulatorExporterCurveFittingModes resolve_curve_fitting_mode(CircuitSimulatorExporterCurveFittingModes mode);
        static std::vector<double> calculate_core_resistance_coefficients(Magnetic magnetic, double temperature = defaults.ambientTemperature, CoreLossTopology topology = CoreLossTopology::RIDLEY);
        std::string export_magnetic_as_symbol(Magnetic magnetic, std::optional<std::string> outputFilename = std::nullopt, std::optional<std::string> filePathOrFile = std::nullopt);
        std::string export_magnetic_as_subcircuit(Magnetic magnetic, double frequency = defaults.measurementFrequency, double temperature = defaults.ambientTemperature, std::optional<std::string> outputFilename = std::nullopt, std::optional<std::string> filePathOrFile = std::nullopt, CircuitSimulatorExporterCurveFittingModes mode=CircuitSimulatorExporterCurveFittingModes::LADDER);
};

class CircuitSimulatorExporterSimbaModel : public CircuitSimulatorExporterModel {
    public:
        std::string programName = "Simba";
        double _scale;
        double _modelSize = 50;
        size_t _precision = 12;

        std::string generate_id();
        std::mt19937 _gen;
        CircuitSimulatorExporterSimbaModel();
        std::string export_magnetic_as_symbol(Magnetic magnetic, std::optional<std::string> filePathOrFile = std::nullopt) {
            return "Not supported";
        }
        std::string export_magnetic_as_subcircuit(Magnetic magnetic, double frequency = defaults.measurementFrequency, double temperature = defaults.ambientTemperature, std::optional<std::string> filePathOrFile = std::nullopt, CircuitSimulatorExporterCurveFittingModes mode=CircuitSimulatorExporterCurveFittingModes::LADDER);
        ordered_json create_device(std::string libraryName, std::vector<int> coordinates, int angle, std::string name);
        ordered_json create_air_gap(std::vector<int> coordinates, double area, double length, int angle, std::string name);
        ordered_json create_core(double initialPermeability, std::vector<int> coordinates, double area, double length, int angle, std::string name);
        ordered_json create_winding(size_t numberTurns, std::vector<int> coordinates, int angle, std::string name);
        ordered_json create_pin(std::vector<int> coordinates, int angle, std::string name);
        ordered_json create_resistor(double resistance, std::vector<int> coordinates, int angle, std::string name);
        ordered_json create_inductor(double inductance, std::vector<int> coordinates, int angle, std::string name);
        ordered_json create_capacitor(double capacitance, std::vector<int> coordinates, int angle, std::string name);
        std::pair<std::vector<ordered_json>, std::vector<ordered_json>> create_ladder(std::vector<double> ladderCoefficients, std::vector<int> coordinates, std::string name);
        std::pair<std::vector<ordered_json>, std::vector<ordered_json>> create_rosano_ladder(std::vector<double> coefficients, std::vector<int> coordinates, std::string name);
        std::pair<std::vector<ordered_json>, std::vector<ordered_json>> create_fracpole_ladder(const FractionalPoleNetwork& net, std::vector<int> coordinates, const std::string& name);
        ordered_json create_magnetic_ground(std::vector<int> coordinates, int angle, std::string name);
        ordered_json create_connector(std::vector<int> startingCoordinates, std::vector<int> endingCoordinates, std::string name);
        ordered_json merge_connectors(ordered_json connectors);
};

class CircuitSimulatorExporterNgspiceModel : public CircuitSimulatorExporterModel {
    public:
        std::string programName = "Ngspice";
        std::string export_magnetic_as_symbol(Magnetic magnetic, std::optional<std::string> filePathOrFile = std::nullopt) {
            return "Not supported";
        }
        std::string export_magnetic_as_subcircuit(Magnetic magnetic, double frequency = defaults.measurementFrequency, double temperature = defaults.ambientTemperature, std::optional<std::string> filePathOrFile = std::nullopt, CircuitSimulatorExporterCurveFittingModes mode=CircuitSimulatorExporterCurveFittingModes::LADDER);
};

class CircuitSimulatorExporterLtspiceModel : public CircuitSimulatorExporterModel {
    public:
        std::string programName = "Ltspice";
        std::string export_magnetic_as_subcircuit(Magnetic magnetic, double frequency = defaults.measurementFrequency, double temperature = defaults.ambientTemperature, std::optional<std::string> filePathOrFile = std::nullopt, CircuitSimulatorExporterCurveFittingModes mode=CircuitSimulatorExporterCurveFittingModes::LADDER);
        std::string export_magnetic_as_symbol(Magnetic magnetic, std::optional<std::string> filePathOrFile = std::nullopt);
};

class CircuitSimulatorExporterPlecsModel : public CircuitSimulatorExporterModel {
    public:
        std::string programName = "PLECS";
        double _modelSize = 400;
        size_t _precision = 6;

        CircuitSimulatorExporterPlecsModel() = default;

        std::string export_magnetic_as_symbol(Magnetic magnetic, std::optional<std::string> filePathOrFile = std::nullopt);
        std::string export_magnetic_as_subcircuit(Magnetic magnetic, double frequency = defaults.measurementFrequency, double temperature = defaults.ambientTemperature, std::optional<std::string> filePathOrFile = std::nullopt, CircuitSimulatorExporterCurveFittingModes mode = CircuitSimulatorExporterCurveFittingModes::LADDER);

    private:
        std::string encode_init_commands(const std::string& plainText);
        std::string emit_header(const std::string& name);
        std::string emit_footer();
        std::string emit_schematic_header();
        std::string emit_schematic_footer();
        std::string emit_magnetic_interface(const std::string& name, const std::string& turnsVar, int polarity, std::vector<int> position, const std::string& direction, bool flipped);
        std::string emit_p_sat(const std::string& name, const std::string& areaVar, const std::string& lengthVar, const std::string& muRVar, const std::string& bSatVar, std::vector<int> position, const std::string& direction, bool flipped);
        std::string emit_p_air(const std::string& name, const std::string& areaVar, const std::string& lengthVar, std::vector<int> position, const std::string& direction, bool flipped);
        std::string emit_ac_voltage_source(const std::string& name, std::vector<int> position);
        std::string emit_resistor(const std::string& name, const std::string& valueVar, std::vector<int> position);
        std::string emit_scope(const std::string& name, std::vector<int> position);
        std::string emit_probe(const std::string& name, std::vector<int> position, const std::vector<std::pair<std::string, std::string>>& probes);
        std::string emit_connection(const std::string& type, const std::string& srcComponent, int srcTerminal, const std::string& dstComponent, int dstTerminal, std::vector<std::vector<int>> points = {});

        // Shared helpers for subcircuit/symbol export
        void build_physical_init(std::ostringstream& init, Core& core,
                                  const std::vector<Winding>& windings,
                                  const std::vector<ColumnElement>& columns,
                                  double temperature);
        void build_magnetic_schematic(std::ostringstream& init, std::ostringstream& schematic,
                                       Core& core,
                                       const std::vector<Winding>& windings,
                                       const std::vector<ColumnElement>& columns,
                                       bool isMultiColumn, int yBase);
        void build_electrical_schematic(std::ostringstream& schematic,
                                         const std::vector<Winding>& windings,
                                         bool isMultiColumn, int yBase);
        std::string assemble_plecs_file(const std::string& name, const std::string& initStr, const std::string& schematicStr);
};

class CircuitSimulationReader {
  public:
    enum class DataType : int {
        TIME,
        VOLTAGE,
        CURRENT,
        MAGNETIZING_CURRENT,
        UNKNOWN
    };

  private:
    class CircuitSimulationSignal {
        public:
            CircuitSimulationSignal() = default;
            virtual ~CircuitSimulationSignal() = default;
            std::string name;
            std::vector<double> data;
            DataType type;
            size_t windingIndex;
            size_t operatingPointIndex;

        public:

    };


    std::vector<CircuitSimulationSignal> _columns;
    std::vector<Waveform> _waveforms;
    CircuitSimulationSignal _time;

    std::optional<size_t> _periodStartIndex = std::nullopt;
    std::optional<size_t> _periodStopIndex = std::nullopt;

    std::vector<std::string> _timeAliases = {"TIME", "Time", "time", "[s]"};
    std::vector<std::string> _magnetizingCurrentAliases = {"Imag", "MAG", "mag", "Im", "magnetizing"};
    std::vector<std::string> _currentAliases = {"Inductor current", "CURRENT", "CURR", "Current", "Curr", "current", "curr", "Ipri", "Isec", "I_", "i_", "I(", "i(", "[A]", "Ip", "Is", "It", "Id"};
    std::vector<std::string> _voltageAliases = {"Inductor voltage", "VOLTAGE", "VOLT", "Voltage", "Volt", "voltage", "volt", "Vpri", "Vsec", "Vout", "V_", "v_", "V(", "v(", "[V]", "Vp", "Vs", "Vt"};
    std::vector<std::string> _currentPrefixes = {"I"};
    std::vector<std::string> _voltagePrefixes = {"V"};

  public:

    CircuitSimulationReader() = default;
    virtual ~CircuitSimulationReader() = default;

    CircuitSimulationReader(std::string filePathOrFile, bool forceFile=false);

    void process_line(std::string line, char separator);
    bool extract_winding_indexes(size_t numberWindings);
    bool extract_column_types(double frequency);
    OperatingPoint extract_operating_point(size_t numberWindings, double frequency, std::optional<std::vector<std::map<std::string, std::string>>> mapColumnNames = std::nullopt, double ambientTemperature = 25);
    std::vector<std::map<std::string, std::string>> extract_map_column_names(size_t numberWindings, double frequency);
    bool assign_column_names(std::vector<std::map<std::string, std::string>> columnNames);
    std::vector<std::string> extract_column_names();
    Waveform extract_waveform(CircuitSimulationSignal signal, double frequency, bool sample=true);
    static CircuitSimulationSignal find_time(std::vector<CircuitSimulationSignal> columns);
    Waveform get_one_period(Waveform waveform, double frequency, bool sample=true, bool alignToZeroCrossing=true);
    static char guess_separator(std::string line);
    static bool can_be_voltage(std::vector<double> data, double limit=0.05);
    static bool can_be_time(std::vector<double> data);
    static bool can_be_current(std::vector<double> data, double limit=0.05);
    std::optional<CircuitSimulationReader::DataType> guess_type_by_name(std::string name);
};

class CircuitSimulatorExporterNl5Model {
    private:
        int _precision = 10;

    public:
        CircuitSimulatorExporterNl5Model() = default;

        std::string export_magnetic_as_subcircuit(
            Magnetic magnetic,
            double frequency,
            double temperature,
            std::optional<std::string> filePathOrFile = std::nullopt,
            CircuitSimulatorExporterCurveFittingModes mode =
                CircuitSimulatorExporterCurveFittingModes::AUTO);
};




} // namespace OpenMagnetics
