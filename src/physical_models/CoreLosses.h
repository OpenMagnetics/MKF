#pragma once
#include "Constants.h"
#include "Models.h"

#include "constructive_models/Core.h"
#include "processors/Inputs.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <streambuf>
#include <vector>
#include "support/Exceptions.h"

using namespace MAS;

namespace OpenMagnetics {

inline std::map<std::string, tk::spline> lossFactorInterps;

// ============================================================================
// Roshen Hysteresis Model - Data Structures
// ============================================================================

/**
 * @brief Cache key for B-H loop calculations
 */
struct BHLoopCacheKey {
    std::string materialName;
    double temperature;
    double magneticFluxDensityPeak;
    
    bool operator<(const BHLoopCacheKey& other) const {
        if (materialName != other.materialName) return materialName < other.materialName;
        if (temperature != other.temperature) return temperature < other.temperature;
        return magneticFluxDensityPeak < other.magneticFluxDensityPeak;
    }
    
    bool operator==(const BHLoopCacheKey& other) const {
        return materialName == other.materialName &&
               temperature == other.temperature &&
               magneticFluxDensityPeak == other.magneticFluxDensityPeak;
    }
};

/**
 * @brief Represents a complete B-H hysteresis loop
 */
struct BHLoop {
    std::vector<double> upperH;  // H values for upper branch
    std::vector<double> upperB;  // B values for upper branch
    std::vector<double> lowerH;  // H values for lower branch
    std::vector<double> lowerB;  // B values for lower branch
    double area;                 // Loop area (hysteresis loss per cycle)
    double coerciveForce;
    double remanence;
    double saturationB;
    double saturationH;
    
    /**
     * @brief Calculate loop area using trapezoidal integration
     */
    double calculate_area() const;
    
    /**
     * @brief Get B value for given H on upper or lower branch
     */
    double get_b_for_h(double H, bool upperBranch) const;
};

/**
 * @brief Parameters for Roshen hysteresis model
 */
struct RoshenParameters {
    double coerciveForce;
    double remanence;
    double saturationMagneticFluxDensity;
    double saturationMagneticFieldStrength;
    double a1;  // Roshen coefficient
    double b1;  // Roshen coefficient
    double b2;  // Roshen coefficient
    
    // Temperature coefficients for correction
    double alpha_Hc = 0.002;  // Hc temperature coefficient (1/°C)
    double alpha_Br = 0.001;  // Br temperature coefficient (1/°C)
    double T_ref = 25.0;      // Reference temperature (°C)
    
    /**
     * @brief Calculate Roshen coefficients from material properties
     */
    static RoshenParameters from_material(const CoreMaterial& material, double temperature);
    
    /**
     * @brief Apply temperature correction to coercive force and remanence
     */
    void apply_temperature_correction(double temperature);
    
    /**
     * @brief Calculate squared error between model and measured data
     */
    double calculate_fitting_error(const std::vector<std::pair<double, double>>& measuredPoints) const;
};

/**
 * @brief Improved Roshen hysteresis model with caching and analytical methods
 * 
 * This class provides:
 * - Cached B-H loop calculations for improved performance
 * - Analytical minor loop scaling (O(N) instead of iterative O(N*k))
 * - Sub-loop handling for arbitrary waveforms
 * - DC bias correction
 * - Air gap shearing support
 */
class RoshenHysteresisModel {
private:
    // Cache for computed loops
    static std::map<BHLoopCacheKey, BHLoop> _loopCache;
    static std::map<std::string, RoshenParameters> _paramsCache;
    
    // Cache configuration
    static constexpr size_t MAX_CACHE_SIZE = 1000;
    
    /**
     * @brief Generate H points from -Hs to +Hs
     */
    std::vector<double> generate_h_points(double Hs, double step) const;
    
    /**
     * @brief Calculate B from H using Roshen equation
     */
    double calculate_b(double H, double coerciveForce, double a, double b, bool upperBranch) const;
    
    /**
     * @brief Calculate B-H curve point using Roshen hyperbolic model
     */
    double bh_curve_half_loop(double H, double Hc, double a, double b) const {
        return (H + Hc) / (a + b * std::fabs(H + Hc));
    }

public:
    RoshenHysteresisModel() = default;
    ~RoshenHysteresisModel() = default;
    
    /**
     * @brief Clear all caches
     */
    static void clear_cache();
    
    /**
     * @brief Get cached parameters for a material, or calculate and cache them
     */
    static RoshenParameters get_cached_parameters(const CoreMaterial& material, double temperature);
    
    /**
     * @brief Calculate major B-H loop (full saturation)
     */
    BHLoop calculate_major_loop(const RoshenParameters& params, double temperature, double stepSize = -1.0);
    
    /**
     * @brief Analytically scale major loop to minor loop at target flux density
     * 
     * This replaces the iterative method with O(N) analytical scaling.
     * 
     * @param majorLoop The full major hysteresis loop
     * @param params Material parameters
     * @param targetBpeak Target peak flux density for minor loop
     * @return Scaled minor loop
     */
    BHLoop scale_to_minor_loop_analytical(const BHLoop& majorLoop, 
                                          const RoshenParameters& params,
                                          double targetBpeak);
    
    /**
     * @brief Get minor loop with caching
     */
    BHLoop get_minor_loop(const CoreMaterial& material, 
                         double temperature,
                         double targetBpeak,
                         bool useCache = true);
    
    /**
     * @brief Calculate hysteresis loss density from B-H loop area
     */
    double calculate_hysteresis_loss_density(const BHLoop& loop, double frequency) const {
        return loop.area * frequency;
    }
    
    /**
     * @brief Detect and calculate sub-loops in an arbitrary waveform
     * 
     * Based on Roshen Eq. (8.21): Total loss = sum of losses from each separated loop
     * 
     * @param H_wavefield Magnetic field strength waveform
     * @param B_wavefield Magnetic flux density waveform
     * @param frequency Operating frequency
     * @return Total hysteresis loss density including sub-loops
     */
    double calculate_loss_with_subloops(const std::vector<double>& H_wavefield,
                                       const std::vector<double>& B_wavefield,
                                       double frequency);
    
    /**
     * @brief Apply DC bias correction to loop
     * 
     * DC bias shifts the loop and changes effective permeability
     */
    BHLoop apply_dc_bias(const BHLoop& originalLoop, double B_dc, const RoshenParameters& params);
    
    /**
     * @brief Apply air gap shearing to loop
     * 
     * Air gaps cause shearing of the B-H loop (tilt)
     */
    BHLoop apply_air_gap_shear(const BHLoop& originalLoop, 
                              double gapLength, 
                              double coreLength,
                              double effectivePermeability);
    
    /**
     * @brief High-level interface matching existing API
     * 
     * Returns upper and lower B waveforms for given H field
     */
    std::pair<std::vector<double>, std::vector<double>> get_hysteresis_loop(
        const CoreMaterial& material,
        double temperature,
        std::optional<double> magneticFluxDensityPeak,
        std::optional<double> magneticFieldStrengthPeak);
    
    /**
     * @brief Get amplitude permeability from minor loop slope
     */
    std::optional<double> get_amplitude_permeability(const CoreMaterial& material,
                                                     double temperature,
                                                     std::optional<double> magneticFluxDensityPeak,
                                                     std::optional<double> magneticFieldStrengthPeak);
    
    /**
     * @brief Fit Roshen parameters to measured B-H data using least squares
     * 
     * This provides better accuracy than analytical calculation from saturation points alone.
     * Uses the existing levmar library for optimization.
     * 
     * @param measuredPoints Vector of (H, B) measured data points
     * @param initialParams Initial parameter estimate
     * @return Optimized parameters
     */
    static RoshenParameters fit_parameters_to_data(
        const std::vector<std::pair<double, double>>& measuredPoints,
        const RoshenParameters& initialParams);
    
    /**
     * @brief Calculate loop area using adaptive Simpson integration
     * 
     * More accurate than simple trapezoidal rule, especially near coercive points
     * where curvature is high.
     * 
     * @param loop B-H loop to integrate
     * @param tolerance Integration tolerance
     * @return Loop area
     */
    double calculate_loop_area_adaptive(const BHLoop& loop, double tolerance = 1e-6) const;
    
    /**
     * @brief Get parameters with temperature correction applied
     * 
     * @param material Core material
     * @param temperature Operating temperature
     * @return Parameters with temperature-corrected Hc and Br
     */
    static RoshenParameters get_temperature_corrected_parameters(
        const CoreMaterial& material, 
        double temperature);
    
    /**
     * @brief Calculate hysteresis losses with all improvements
     * 
     * Combines temperature correction, dynamic fitting (if data available),
     * and adaptive integration for maximum accuracy.
     * 
     * @param core Core geometry
     * @param excitation Operating point excitation
     * @param temperature Core temperature
     * @return Hysteresis loss density (W/m³)
     */
    double calculate_hysteresis_losses_improved(
        Core core,
        OperatingPointExcitation excitation,
        double temperature);
};

// ============================================================================
// Core Losses Models
// ============================================================================

class CoreLossesModel {
  private:
  protected:
    double _get_frequency_from_core_losses(Core core,
                                           SignalDescriptor magneticFluxDensity,
                                           double temperature,
                                           double coreLosses);
    SignalDescriptor _get_magnetic_flux_density_from_core_losses(Core core,
                                                                        double frequency,
                                                                        double temperature,
                                                                        double coreLosses);

  public:
    std::vector<double> _hysteresisMajorLoopTop;
    std::vector<double> _hysteresisMajorLoopBottom;
    std::vector<double> _hysteresisMajorH;
    std::vector<double> _hysteresisMinorLoopTop;
    std::vector<double> _hysteresisMinorLoopBottom;
    std::vector<double> _hysteresisMinorH;
    SteinmetzCoreLossesMethodRangeDatum _steinmetzDatum;
    bool _steinmetzDatumSet = false;
    CoreLossesModel() = default;
    virtual ~CoreLossesModel() = default;
    virtual CoreLossesOutput get_core_losses(Core core,
                                             OperatingPointExcitation excitation,
                                             double temperature) = 0;
    virtual double get_core_volumetric_losses(CoreMaterial coreMaterial,
                                             OperatingPointExcitation excitation,
                                             double temperature) = 0;
    virtual double get_core_mass_losses(CoreMaterial coreMaterial,
                                             OperatingPointExcitation excitation,
                                             double temperature) = 0;
    virtual double get_frequency_from_core_losses(Core core,
                                                  SignalDescriptor magneticFluxDensity,
                                                  double temperature,
                                                  double coreLosses) = 0;
    virtual SignalDescriptor get_magnetic_flux_density_from_core_losses(Core core,
                                                                                double frequency,
                                                                                double temperature,
                                                                                double coreLosses) = 0;
    static std::shared_ptr<CoreLossesModel> factory(CoreLossesModels modelName);
    static std::shared_ptr<CoreLossesModel> factory(std::map<std::string, std::string> models);
    static std::shared_ptr<CoreLossesModel> factory(json models);
    static SteinmetzCoreLossesMethodRangeDatum get_steinmetz_coefficients(
        CoreMaterialDataOrNameUnion material,
        double frequency);
    bool is_steinmetz_datum_loaded() { return _steinmetzDatumSet; }
    static CoreLossesMethodData get_method_data(CoreMaterial materialData, std::string method);
    static std::vector<VolumetricLossesPoint> get_volumetric_losses_data(CoreMaterial materialData);
    SteinmetzCoreLossesMethodRangeDatum get_steinmetz_datum() { return _steinmetzDatum; }
    void set_steinmetz_datum(SteinmetzCoreLossesMethodRangeDatum steinmetzDatum) {
        _steinmetzDatumSet = true;
        _steinmetzDatum = steinmetzDatum;
    }

    double get_core_losses_series_resistance(Core core,
                                             double frequency,
                                             double temperature,
                                             double magnetizingInductance);

    bool usesVolumetricLosses(CoreMaterial material) {
        if (material.get_volumetric_losses().size() != 0) {
            return true;
        }
        else {
            if (!material.get_mass_losses()) {
                throw std::runtime_error("Material is missing voluemtric and massic losses");
            }
            return false;
        }
    }

    double get_magnetic_flux_density_from_volumetric_losses(SteinmetzCoreLossesMethodRangeDatum steinmetzDatum, double volumetricLosses, double frequency, double temperature) {
        double temperatureTerm = 1;
        double k = steinmetzDatum.get_k();
        double alpha = steinmetzDatum.get_alpha();
        double beta = steinmetzDatum.get_beta();
        if (steinmetzDatum.get_ct0() && steinmetzDatum.get_ct1() && steinmetzDatum.get_ct2()) {
            double ct0 = steinmetzDatum.get_ct0().value();
            double ct1 = steinmetzDatum.get_ct1().value();
            double ct2 = steinmetzDatum.get_ct2().value();
            temperatureTerm = ct2 * pow(temperature, 2) - ct1 * temperature + ct0;
        }
        return pow(volumetricLosses / k / pow(frequency, alpha) / temperatureTerm, 1 / beta);
    }

    static double apply_temperature_coefficients(double volumetricLosses,
                                                 SteinmetzCoreLossesMethodRangeDatum steinmetzDatum,
                                                 double temperature) {
        double volumetricLossesWithTemperature = volumetricLosses;
        if (steinmetzDatum.get_ct0() && steinmetzDatum.get_ct1() && steinmetzDatum.get_ct2()) {
            double ct0 = steinmetzDatum.get_ct0().value();
            double ct1 = steinmetzDatum.get_ct1().value();
            double ct2 = steinmetzDatum.get_ct2().value();
            double scale = (ct2 * pow(temperature, 2) - ct1 * temperature + ct0);
            if (scale > 0) {
                volumetricLossesWithTemperature *= scale;
            }
        }
        return volumetricLossesWithTemperature;
    }

    static std::vector<std::string> get_methods_string(CoreMaterialDataOrNameUnion material) {
        std::vector<std::string> methodsString;
        auto methods = get_methods(material);
        for (auto method : methods) {
            json methodString;
            to_json(methodString, method);
            methodsString.push_back(methodString);
        }
        return methodsString;
    }
    static std::vector<CoreLossesModels> get_methods(CoreMaterialDataOrNameUnion material);

    static std::map<std::string, std::string> get_models_information() {
        std::map<std::string, std::string> information;
        information["Steinmetz"] = R"(Based on "On the law of hysteresis" by Charles Proteus Steinmetz)";
        information["iGSE"] =
            R"(Based on "Accurate Prediction of Ferrite Core Loss with Nonsinusoidal Waveforms Using Only Steinmetz Parameters" by Charles R. Sullivan)";
        information["Barg"] =
            R"(Based on "Core Loss Calculation of Symmetric Trapezoidal Magnetic Flux Density Waveform" by Sobhi Barg)";
        information["Roshen"] =
            R"(Based on "Ferrite Core Loss for Power Magnetic Components Design" and "A Practical, Accurate and Very General Core Loss Model for Nonsinusoidal Waveforms" by Waseem Roshen)";
        information["Albach"] =
            R"(Based on "Calculating Core Losses in Transformers for Arbitrary Magnetizing Currents A Comparison of Different Approaches" by Manfred Albach)";
        information["NSE"] =
            R"(Based on "Measurement and Loss Model of Ferrites with Non-sinusoidal Waveforms" by Alex Van den Bossche)";
        information["MSE"] =
            R"(Based on "Calculation of Losses in Ferro- and Ferrimagnetic Materials Based on the Modified Steinmetz Equation" by Jürgen Reinert)";
        return information;
    }
    static std::map<std::string, double> get_models_errors() {
        // This are taken manually from running the tests. Yes, a pain in the ass
        // TODO: Automate it somehow
        std::map<std::string, double> errors;
        errors["Steinmetz"] = 0.39353;
        errors["iGSE"] = 0.358237;
        errors["Barg"] = 0.374326;
        errors["Roshen"] = 0.487881;
        errors["Albach"] = 0.357267;
        errors["NSE"] = 0.358237;
        errors["MSE"] = 0.357267;
        return errors;
    }
    static std::map<std::string, std::string> get_models_external_links() {
        std::map<std::string, std::string> external_links;
        external_links["Steinmetz"] = "https://ieeexplore.ieee.org/document/1457110";
        external_links["iGSE"] = "http://inductor.thayerschool.org/papers/IGSE.pdf";
        external_links["Barg"] = "https://miun.diva-portal.org/smash/get/diva2:1622559/FULLTEXT01.pdf";
        external_links["Roshen"] = "https://ieeexplore.ieee.org/document/4052433";
        external_links["Albach"] = "https://ieeexplore.ieee.org/iel3/3925/11364/00548774.pdf";
        external_links["NSE"] = "http://web.eecs.utk.edu/~dcostine/ECE482/Spring2015/materials/magnetics/NSE.pdf";
        external_links["MSE"] = "https://www.mikrocontroller.net/attachment/129490/Modified_Steinmetz.pdf";
        return external_links;
    }
    static std::map<std::string, std::string> get_models_internal_links() {
        std::map<std::string, std::string> internal_links;
        internal_links["Steinmetz"] = "";
        internal_links["iGSE"] = "";
        internal_links["Barg"] = "";
        internal_links["Roshen"] = "/musings/4_roshen_method_core_losses";
        internal_links["Albach"] = "";
        internal_links["NSE"] = "";
        internal_links["MSE"] = "";
        return internal_links;
    }
};

// Based on On the law of hysteresis by Charles Proteus Steinmetz
// https://sci-hub.st/10.1109/proc.1984.12842
class CoreLossesSteinmetzModel : public CoreLossesModel {
  private:
    std::string _modelName = "Steinmetz";
  public:
    CoreLossesOutput get_core_losses(Core core,
                                     OperatingPointExcitation excitation,
                                     double temperature);
    double get_core_volumetric_losses(CoreMaterial coreMaterial,
                                      OperatingPointExcitation excitation,
                                      double temperature);
    double get_core_mass_losses(CoreMaterial coreMaterial,
                                      OperatingPointExcitation excitation,
                                      double temperature) {
        throw std::runtime_error("Mass losses is only valid for Proprietary models from Magnetec");
    }
    double get_frequency_from_core_losses(Core core,
                                          SignalDescriptor magneticFluxDensity,
                                          double temperature,
                                          double coreLosses);

    static std::pair<std::vector<SteinmetzCoreLossesMethodRangeDatum>, std::vector<double>> calculate_steinmetz_coefficients(std::vector<VolumetricLossesPoint> volumetricLosses, std::vector<std::pair<double, double>> ranges);

    SignalDescriptor get_magnetic_flux_density_from_core_losses(Core core,
                                                                        double frequency,
                                                                        double temperature,
                                                                        double coreLosses);
};

// Based on Accurate Prediction of Ferrite Core Loss with Nonsinusoidal Waveforms Using Only Steinmetz Parameters by
// Charles R. Sullivan http://inductor.thayerschool.org/papers/IGSE.pdf
class CoreLossesIGSEModel : public CoreLossesSteinmetzModel {
  private:
    std::string _modelName = "iGSE";
  public:
    double get_core_volumetric_losses(CoreMaterial coreMaterial,
                                     OperatingPointExcitation excitation,
                                     double temperature);
    double get_core_mass_losses(CoreMaterial coreMaterial,
                                      OperatingPointExcitation excitation,
                                      double temperature) {
        throw std::runtime_error("Mass losses is only valid for Proprietary models from Magnetec");
    }
    double get_frequency_from_core_losses(Core core,
                                          SignalDescriptor magneticFluxDensity,
                                          double temperature,
                                          double coreLosses) {
        return _get_frequency_from_core_losses(core, magneticFluxDensity, temperature, coreLosses);
    }
    SignalDescriptor get_magnetic_flux_density_from_core_losses(Core core,
                                                                        double frequency,
                                                                        double temperature,
                                                                        double coreLosses) {
        return _get_magnetic_flux_density_from_core_losses(core, frequency, temperature, coreLosses);
    }
    double get_ki(SteinmetzCoreLossesMethodRangeDatum steinmetzDatum);
};

// Based on Core Loss Calculation of Symmetric Trapezoidal Magnetic Flux Density Waveform by Sobhi Barg
// https://miun.diva-portal.org/smash/get/diva2:1622559/FULLTEXT01.pdf
class CoreLossesBargModel : public CoreLossesSteinmetzModel {
  private:
    std::string _modelName = "Barg";
  public:
    double get_core_volumetric_losses(CoreMaterial coreMaterial,
                                     OperatingPointExcitation excitation,
                                     double temperature);
    double get_core_mass_losses(CoreMaterial coreMaterial,
                                      OperatingPointExcitation excitation,
                                      double temperature) {
        throw std::runtime_error("Mass losses is only valid for Proprietary models from Magnetec");
    }
    double get_frequency_from_core_losses(Core core,
                                          SignalDescriptor magneticFluxDensity,
                                          double temperature,
                                          double coreLosses) {
        return _get_frequency_from_core_losses(core, magneticFluxDensity, temperature, coreLosses);
    }
    SignalDescriptor get_magnetic_flux_density_from_core_losses(Core core,
                                                                        double frequency,
                                                                        double temperature,
                                                                        double coreLosses) {
        return _get_magnetic_flux_density_from_core_losses(core, frequency, temperature, coreLosses);
    }
};

// Based on Ferrite Core Loss for Power Magnetic Components Design and
// A Practical, Accurate and Very General Core Loss Model for Nonsinusoidal Waveforms by Waseem Roshen
// https://sci-hub.st/10.1109/20.278656
// https://sci-hub.st/10.1109/TPEL.2006.886608
class CoreLossesRoshenModel : public CoreLossesModel {
  private:
    std::string _modelName = "Roshen";
  public:
    CoreLossesOutput get_core_losses(Core core,
                                     OperatingPointExcitation excitation,
                                     double temperature);
    double get_core_volumetric_losses(CoreMaterial coreMaterial,
                                     OperatingPointExcitation excitation,
                                     double temperature);
    double get_core_mass_losses(CoreMaterial coreMaterial,
                                      OperatingPointExcitation excitation,
                                      double temperature) {
        throw std::runtime_error("Mass losses is only valid for Proprietary models from Magnetec");
    }
    double get_frequency_from_core_losses(Core core,
                                          SignalDescriptor magneticFluxDensity,
                                          double temperature,
                                          double coreLosses) {
        return _get_frequency_from_core_losses(core, magneticFluxDensity, temperature, coreLosses);
    }

    SignalDescriptor get_magnetic_flux_density_from_core_losses(Core core,
                                                                        double frequency,
                                                                        double temperature,
                                                                        double coreLosses) {
        return _get_magnetic_flux_density_from_core_losses(core, frequency, temperature, coreLosses);
    }
    std::pair<std::vector<double>, std::vector<double>> get_bh_loop(std::map<std::string, double> parameters, OperatingPointExcitation excitation);
    double get_hysteresis_losses_density(std::map<std::string, double> parameters, OperatingPointExcitation excitation);
    double get_eddy_current_losses_density(Core core, OperatingPointExcitation excitation, double resistivity);
    double get_excess_eddy_current_losses_density(OperatingPointExcitation excitation,
                                                  double resistivity,
                                                  double alphaTimesN0);
    std::map<std::string, double> get_roshen_parameters(Core core,
                                                        OperatingPointExcitation excitation,
                                                        double temperature);
};

// Based on Calculating Core Losses in Transformers for Arbitrary Magnetizing Currents A Comparison of Different
// Approaches by Manfred Albach https://sci-hub.st/10.1109/PESC.1996.548774
class CoreLossesAlbachModel : public CoreLossesSteinmetzModel {
  private:
    std::string _modelName = "Albach";
  public:
    double get_core_volumetric_losses(CoreMaterial coreMaterial,
                                     OperatingPointExcitation excitation,
                                     double temperature);
    double get_core_mass_losses(CoreMaterial coreMaterial,
                                      OperatingPointExcitation excitation,
                                      double temperature) {
        throw std::runtime_error("Mass losses is only valid for Proprietary models from Magnetec");
    }
    double get_frequency_from_core_losses(Core core,
                                          SignalDescriptor magneticFluxDensity,
                                          double temperature,
                                          double coreLosses) {
        return _get_frequency_from_core_losses(core, magneticFluxDensity, temperature, coreLosses);
    }
    SignalDescriptor get_magnetic_flux_density_from_core_losses(Core core,
                                                                        double frequency,
                                                                        double temperature,
                                                                        double coreLosses) {
        return _get_magnetic_flux_density_from_core_losses(core, frequency, temperature, coreLosses);
    }
};

// Based on Measurement and Loss Model of Ferrites with Non-sinusoidal Waveforms by Alex Van den Bossche
// http://web.eecs.utk.edu/~dcostine/ECE482/Spring2015/materials/magnetics/NSE.pdf
class CoreLossesNSEModel : public CoreLossesSteinmetzModel {
  private:
    std::string _modelName = "NSE";
  public:
    double get_core_volumetric_losses(CoreMaterial coreMaterial,
                                     OperatingPointExcitation excitation,
                                     double temperature);
    double get_core_mass_losses(CoreMaterial coreMaterial,
                                      OperatingPointExcitation excitation,
                                      double temperature) {
        throw std::runtime_error("Mass losses is only valid for Proprietary models from Magnetec");
    }
    double get_frequency_from_core_losses(Core core,
                                          SignalDescriptor magneticFluxDensity,
                                          double temperature,
                                          double coreLosses) {
        return _get_frequency_from_core_losses(core, magneticFluxDensity, temperature, coreLosses);
    }
    SignalDescriptor get_magnetic_flux_density_from_core_losses(Core core,
                                                                        double frequency,
                                                                        double temperature,
                                                                        double coreLosses) {
        return _get_magnetic_flux_density_from_core_losses(core, frequency, temperature, coreLosses);
    }
    double get_kn(SteinmetzCoreLossesMethodRangeDatum steinmetzDatum);
};

// Based on Calculation of Losses in Ferro- and Ferrimagnetic Materials Based on the Modified Steinmetz Equation by
// Jürgen Reinert https://sci-hub.st/10.1109/28.936396
class CoreLossesMSEModel : public CoreLossesSteinmetzModel {
  private:
    std::string _modelName = "MSE";
  public:
    double get_core_volumetric_losses(CoreMaterial coreMaterial,
                                     OperatingPointExcitation excitation,
                                     double temperature);
    double get_core_mass_losses(CoreMaterial coreMaterial,
                                      OperatingPointExcitation excitation,
                                      double temperature) {
        throw std::runtime_error("Mass losses is only valid for Proprietary models from Magnetec");
    }
    double get_frequency_from_core_losses(Core core,
                                          SignalDescriptor magneticFluxDensity,
                                          double temperature,
                                          double coreLosses) {
        return _get_frequency_from_core_losses(core, magneticFluxDensity, temperature, coreLosses);
    }
    SignalDescriptor get_magnetic_flux_density_from_core_losses(Core core,
                                                                        double frequency,
                                                                        double temperature,
                                                                        double coreLosses) {
        return _get_magnetic_flux_density_from_core_losses(core, frequency, temperature, coreLosses);
    }
};

// Based on Improved Calculation of Core Loss With Nonsinusoidal Waveforms by Charles R. Sullivan
// http://inductor.thayerschool.org/papers/gse.pdf
class CoreLossesGSEModel : public CoreLossesSteinmetzModel {
  private:
    std::string _modelName = "GSE";
  public:
    double get_core_volumetric_losses(CoreMaterial coreMaterial,
                                     OperatingPointExcitation excitation,
                                     double temperature);
    double get_core_mass_losses(CoreMaterial coreMaterial,
                                      OperatingPointExcitation excitation,
                                      double temperature) {
        throw std::runtime_error("Mass losses is only valid for Proprietary models from Magnetec");
    }
    double get_frequency_from_core_losses(Core core,
                                          SignalDescriptor magneticFluxDensity,
                                          double temperature,
                                          double coreLosses) {
        return _get_frequency_from_core_losses(core, magneticFluxDensity, temperature, coreLosses);
    }
    SignalDescriptor get_magnetic_flux_density_from_core_losses(Core core,
                                                                        double frequency,
                                                                        double temperature,
                                                                        double coreLosses) {
        return _get_magnetic_flux_density_from_core_losses(core, frequency, temperature, coreLosses);
    }
};

// Based on the formula provided by the manufacturer
class CoreLossesProprietaryModel : public CoreLossesSteinmetzModel {
  private:
    std::string _modelName = "Proprietary";
  public:
    double get_core_volumetric_losses(CoreMaterial coreMaterial,
                                     OperatingPointExcitation excitation,
                                     double temperature);
    double get_core_mass_losses(CoreMaterial coreMaterial,
                                     OperatingPointExcitation excitation,
                                     double temperature);
    double get_frequency_from_core_losses(Core core,
                                          SignalDescriptor magneticFluxDensity,
                                          [[maybe_unused]]double temperature,
                                          double coreLosses);
    SignalDescriptor get_magnetic_flux_density_from_core_losses(Core core,
                                                                        double frequency,
                                                                        double temperature,
                                                                        double coreLosses) {
        return _get_magnetic_flux_density_from_core_losses(core, frequency, temperature, coreLosses);
    }
    static std::map<std::string, std::string> get_core_volumetric_losses_equations(CoreMaterial coreMaterial);
    static std::map<std::string, std::string> get_core_volumetric_losses_equations(CoreLossesMethodData coreLossesMethodData);

};

// Based on the formula provided by the manufacturer
class CoreLossesLossFactorModel : public CoreLossesModel {
  private:
    std::string _modelName = "Loss Factor";
  public:
    CoreLossesOutput get_core_losses(Core core,
                                     OperatingPointExcitation excitation,
                                     double temperature);
    double get_core_volumetric_losses(CoreMaterial coreMaterial,
                                     OperatingPointExcitation excitation,
                                     double temperature) {
        return get_core_volumetric_losses(coreMaterial, excitation, temperature, 1);
    }
    double get_core_mass_losses(CoreMaterial coreMaterial,
                                      OperatingPointExcitation excitation,
                                      double temperature) {
        throw std::runtime_error("Mass losses is only valid for Proprietary models from Magnetec");
    }

    double get_core_losses_series_resistance(CoreMaterial coreMaterial,
                                             double frequency,
                                             double temperature,
                                             double magnetizingInductance);

    double get_core_volumetric_losses(CoreMaterial coreMaterial,
                                     OperatingPointExcitation excitation,
                                     double temperature,
                                     double magnetizingInductance);

    double get_frequency_from_core_losses(Core core,
                                          SignalDescriptor magneticFluxDensity,
                                          double temperature,
                                          double coreLosses) {
        return _get_frequency_from_core_losses(core, magneticFluxDensity, temperature, coreLosses);
    }

    SignalDescriptor get_magnetic_flux_density_from_core_losses(Core core,
                                                                        double frequency,
                                                                        double temperature,
                                                                        double coreLosses) {
        return _get_magnetic_flux_density_from_core_losses(core, frequency, temperature, coreLosses);
    }

    static double calculate_magnetizing_inductance_from_excitation(Core core, OperatingPointExcitation excitation, double temperature);

};

class CoreLosses {
    private:
        std::vector<std::pair<CoreLossesModels, std::shared_ptr<CoreLossesModel>>> _coreLossesModels;
    public:

    CoreLosses() {
        for (auto modelName : settings.get_core_losses_model_names()) {
            _coreLossesModels.push_back(std::pair<CoreLossesModels, std::shared_ptr<CoreLossesModel>>{modelName, CoreLossesModel::factory(modelName)});
        }
    }
    virtual ~CoreLosses() = default;

    void set_core_losses_model_name(CoreLossesModels model) {
        settings.set_core_losses_preferred_model_name(model);
        _coreLossesModels.clear();
        for (auto modelName : settings.get_core_losses_model_names()) {
            _coreLossesModels.push_back(std::pair<CoreLossesModels, std::shared_ptr<CoreLossesModel>>{modelName, CoreLossesModel::factory(modelName)});
        }
    }

    CoreLossesOutput calculate_core_losses(Core core, OperatingPointExcitation excitation, double temperature);
    std::shared_ptr<CoreLossesModel> get_core_losses_model(std::string materialName);
    double get_core_volumetric_losses(CoreMaterial coreMaterial, OperatingPointExcitation excitation, double temperature);
    double get_core_losses_series_resistance(Core core, double frequency, double temperature, double magnetizingInductance);
};

} // namespace OpenMagnetics