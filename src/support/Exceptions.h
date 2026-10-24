/**
 * @file Exceptions.h
 * @brief Exception hierarchy for OpenMagnetics library
 * 
 * This file defines a comprehensive exception hierarchy for the OpenMagnetics
 * library, providing structured error handling with descriptive error codes
 * and messages.
 */

#pragma once

#include <exception>
#include <string>
#include <sstream>
#include <optional>
#include <utility>

namespace OpenMagnetics {

/**
 * @brief Error codes for OpenMagnetics exceptions
 */
enum class ErrorCode : int {
    // General errors (1-99)
    UNKNOWN_ERROR = 1,
    NOT_IMPLEMENTED = 2,
    INVALID_ARGUMENT = 3,
    
    // Core related errors (100-199)
    CORE_NOT_PROCESSED = 100,
    CORE_MATERIAL_NOT_FOUND = 101,
    CORE_SHAPE_NOT_FOUND = 102,
    CORE_INVALID_GAPPING = 103,
    CORE_INVALID_GEOMETRY = 104,
    INVALID_CORE_DATA = 105,
    INVALID_CORE_MATERIAL_DATA = 106,
    
    // Coil related errors (200-299)
    COIL_NOT_PROCESSED = 200,
    COIL_WINDING_ERROR = 201,
    COIL_INVALID_TURNS = 202,
    COIL_WIRE_NOT_FOUND = 203,
    COIL_SHORTED_TURNS = 204,
    
    // Material related errors (300-399)
    MATERIAL_NOT_FOUND = 300,
    MATERIAL_DATA_MISSING = 301,
    MATERIAL_INVALID_PROPERTY = 302,
    MATERIAL_FREQUENCY_OUT_OF_SPAN = 303,  // f outside every fitted loss range (ABT #1456)
    MATERIAL_ABOVE_CURIE_TEMPERATURE = 304,  // loss asked of a material hotter than its Curie point (ABT #1485)
    MATERIAL_REQUESTED_LOSS_MODEL_NOT_AVAILABLE = 305,  // an explicitly requested core-loss model the material cannot run (ABT #1497)
    MATERIAL_TEMPERATURE_OUT_OF_RANGE = 306,  // mu_i(T) asked below Curie but outside the tabulated temperatures (ABT #1739)
    MATERIAL_LOSS_FACTOR_ABOVE_SMALL_SIGNAL_FLUX = 307,  // a small-signal loss factor asked at a larger flux (ABT #2006)
    
    // Calculation errors (400-499)
    CALCULATION_NAN_RESULT = 400,
    CALCULATION_DIVERGED = 401,
    CALCULATION_INVALID_INPUT = 402,
    CALCULATION_ERROR = 403,
    CALCULATION_DIVERGENCE = 404,
    
    // Input/Output errors (500-599)
    IO_FILE_NOT_FOUND = 500,
    IO_PARSE_ERROR = 501,
    IO_SCHEMA_VALIDATION_FAILED = 502,
    MISSING_DATA = 503,
    INVALID_INPUT = 504,
    
    // Gap related errors (600-699)
    GAP_NOT_PROCESSED = 600,
    GAP_INVALID_DIMENSIONS = 601,
    GAP_SHAPE_NOT_SET = 602,
    
    // Model errors (700-799)
    MODEL_NOT_AVAILABLE = 700,
    MODEL_INVALID_PARAMETERS = 701,
    
    // Wire related errors (800-899)
    INVALID_WIRE_DATA = 800,
    
    // Bobbin related errors (900-999)
    INVALID_BOBBIN_DATA = 900,
    
    // Coil configuration errors (1000-1099)
    INVALID_COIL_CONFIGURATION = 1000,
    
    // Insulation related errors (1100-1199)
    INVALID_INSULATION_DATA = 1100,
    
    // Design requirements errors (1200-1299)
    INVALID_DESIGN_REQUIREMENTS = 1200,
    
    // Calculation result errors (1300-1399)
    CALCULATION_INVALID_RESULT = 1300,
    CALCULATION_TIMEOUT = 1301,
    INVALID_LOSS_DATA = 1302,
    CALCULATION_CORE_SATURATES = 1303
};

/**
 * @brief Convert error code to string
 */
inline std::string to_string(ErrorCode code) {
    switch (code) {
        case ErrorCode::UNKNOWN_ERROR: return "UNKNOWN_ERROR";
        case ErrorCode::NOT_IMPLEMENTED: return "NOT_IMPLEMENTED";
        case ErrorCode::INVALID_ARGUMENT: return "INVALID_ARGUMENT";
        case ErrorCode::CORE_NOT_PROCESSED: return "CORE_NOT_PROCESSED";
        case ErrorCode::CORE_MATERIAL_NOT_FOUND: return "CORE_MATERIAL_NOT_FOUND";
        case ErrorCode::CORE_SHAPE_NOT_FOUND: return "CORE_SHAPE_NOT_FOUND";
        case ErrorCode::CORE_INVALID_GAPPING: return "CORE_INVALID_GAPPING";
        case ErrorCode::CORE_INVALID_GEOMETRY: return "CORE_INVALID_GEOMETRY";
        case ErrorCode::INVALID_CORE_DATA: return "INVALID_CORE_DATA";
        case ErrorCode::INVALID_CORE_MATERIAL_DATA: return "INVALID_CORE_MATERIAL_DATA";
        case ErrorCode::COIL_NOT_PROCESSED: return "COIL_NOT_PROCESSED";
        case ErrorCode::COIL_WINDING_ERROR: return "COIL_WINDING_ERROR";
        case ErrorCode::COIL_INVALID_TURNS: return "COIL_INVALID_TURNS";
        case ErrorCode::COIL_WIRE_NOT_FOUND: return "COIL_WIRE_NOT_FOUND";
        case ErrorCode::COIL_SHORTED_TURNS: return "COIL_SHORTED_TURNS";
        case ErrorCode::MATERIAL_NOT_FOUND: return "MATERIAL_NOT_FOUND";
        case ErrorCode::MATERIAL_DATA_MISSING: return "MATERIAL_DATA_MISSING";
        case ErrorCode::MATERIAL_INVALID_PROPERTY: return "MATERIAL_INVALID_PROPERTY";
        case ErrorCode::MATERIAL_FREQUENCY_OUT_OF_SPAN: return "MATERIAL_FREQUENCY_OUT_OF_SPAN";
        case ErrorCode::MATERIAL_ABOVE_CURIE_TEMPERATURE: return "MATERIAL_ABOVE_CURIE_TEMPERATURE";
        case ErrorCode::MATERIAL_REQUESTED_LOSS_MODEL_NOT_AVAILABLE: return "MATERIAL_REQUESTED_LOSS_MODEL_NOT_AVAILABLE";
        case ErrorCode::MATERIAL_LOSS_FACTOR_ABOVE_SMALL_SIGNAL_FLUX: return "MATERIAL_LOSS_FACTOR_ABOVE_SMALL_SIGNAL_FLUX";
        case ErrorCode::MATERIAL_TEMPERATURE_OUT_OF_RANGE: return "MATERIAL_TEMPERATURE_OUT_OF_RANGE";
        case ErrorCode::CALCULATION_NAN_RESULT: return "CALCULATION_NAN_RESULT";
        case ErrorCode::CALCULATION_DIVERGED: return "CALCULATION_DIVERGED";
        case ErrorCode::CALCULATION_INVALID_INPUT: return "CALCULATION_INVALID_INPUT";
        case ErrorCode::CALCULATION_ERROR: return "CALCULATION_ERROR";
        case ErrorCode::CALCULATION_DIVERGENCE: return "CALCULATION_DIVERGENCE";
        case ErrorCode::IO_FILE_NOT_FOUND: return "IO_FILE_NOT_FOUND";
        case ErrorCode::IO_PARSE_ERROR: return "IO_PARSE_ERROR";
        case ErrorCode::IO_SCHEMA_VALIDATION_FAILED: return "IO_SCHEMA_VALIDATION_FAILED";
        case ErrorCode::MISSING_DATA: return "MISSING_DATA";
        case ErrorCode::INVALID_INPUT: return "INVALID_INPUT";
        case ErrorCode::GAP_NOT_PROCESSED: return "GAP_NOT_PROCESSED";
        case ErrorCode::GAP_INVALID_DIMENSIONS: return "GAP_INVALID_DIMENSIONS";
        case ErrorCode::GAP_SHAPE_NOT_SET: return "GAP_SHAPE_NOT_SET";
        case ErrorCode::MODEL_NOT_AVAILABLE: return "MODEL_NOT_AVAILABLE";
        case ErrorCode::MODEL_INVALID_PARAMETERS: return "MODEL_INVALID_PARAMETERS";
        case ErrorCode::INVALID_WIRE_DATA: return "INVALID_WIRE_DATA";
        case ErrorCode::INVALID_BOBBIN_DATA: return "INVALID_BOBBIN_DATA";
        case ErrorCode::INVALID_COIL_CONFIGURATION: return "INVALID_COIL_CONFIGURATION";
        case ErrorCode::INVALID_INSULATION_DATA: return "INVALID_INSULATION_DATA";
        case ErrorCode::INVALID_DESIGN_REQUIREMENTS: return "INVALID_DESIGN_REQUIREMENTS";
        case ErrorCode::CALCULATION_INVALID_RESULT: return "CALCULATION_INVALID_RESULT";
        case ErrorCode::CALCULATION_TIMEOUT: return "CALCULATION_TIMEOUT";
        case ErrorCode::INVALID_LOSS_DATA: return "INVALID_LOSS_DATA";
        case ErrorCode::CALCULATION_CORE_SATURATES: return "CALCULATION_CORE_SATURATES";
        default: return "UNKNOWN";
    }
}

/**
 * @brief Base exception class for all OpenMagnetics exceptions
 */
class OpenMagneticsException : public std::exception {
private:
    ErrorCode _code;
    std::string _message;
    std::string _fullMessage;
    std::string _context;
    
public:
    OpenMagneticsException(ErrorCode code, const std::string& message, const std::string& context = "")
        : _code(code), _message(message), _context(context) {
        std::ostringstream oss;
        oss << "[" << to_string(code) << "] " << message;
        if (!context.empty()) {
            oss << " (Context: " << context << ")";
        }
        _fullMessage = oss.str();
    }
    
    virtual ~OpenMagneticsException() noexcept = default;
    
    const char* what() const noexcept override {
        return _fullMessage.c_str();
    }
    
    ErrorCode code() const noexcept { return _code; }
    const std::string& message() const noexcept { return _message; }
    const std::string& context() const noexcept { return _context; }
};

// ============================================================================
// Core Exceptions
// ============================================================================

class CoreException : public OpenMagneticsException {
public:
    CoreException(ErrorCode code, const std::string& message, const std::string& context = "")
        : OpenMagneticsException(code, message, context) {}
};

class CoreNotProcessedException : public CoreException {
public:
    explicit CoreNotProcessedException(const std::string& context = "")
        : CoreException(ErrorCode::CORE_NOT_PROCESSED, 
                       "Core has not been processed. Call process_data() first.", context) {}
};

class CoreMaterialNotFoundException : public CoreException {
public:
    explicit CoreMaterialNotFoundException(const std::string& materialName)
        : CoreException(ErrorCode::CORE_MATERIAL_NOT_FOUND, 
                       "Core material not found: " + materialName, materialName) {}
};

class CoreShapeNotFoundException : public CoreException {
public:
    explicit CoreShapeNotFoundException(const std::string& shapeName)
        : CoreException(ErrorCode::CORE_SHAPE_NOT_FOUND, 
                       "Core shape not found: " + shapeName, shapeName) {}
};

class InvalidGappingException : public CoreException {
public:
    explicit InvalidGappingException(const std::string& message)
        : CoreException(ErrorCode::CORE_INVALID_GAPPING, message) {}
};

// ============================================================================
// Coil Exceptions
// ============================================================================

class CoilException : public OpenMagneticsException {
public:
    CoilException(ErrorCode code, const std::string& message, const std::string& context = "")
        : OpenMagneticsException(code, message, context) {}
};

class CoilNotProcessedException : public CoilException {
public:
    explicit CoilNotProcessedException(const std::string& context = "")
        : CoilException(ErrorCode::COIL_NOT_PROCESSED, 
                       "Coil has not been processed.", context) {}
};

class WireNotFoundException : public CoilException {
public:
    explicit WireNotFoundException(const std::string& wireName)
        : CoilException(ErrorCode::COIL_WIRE_NOT_FOUND,
                       "Wire not found: " + wireName, wireName) {}
};

// ABT #406: two conductors whose COPPER surfaces touch are shorted turns — an electrical fault in
// the described part, not a numerical problem. It surfaced as a divergent turn-to-turn capacitance
// (ABT #395), but the infinity is merely the arithmetic that the fault happens to produce, so name
// the fault instead: a reader told "capacitance is inf" goes looking at the capacitance model,
// while a reader told "these turns are in contact" goes and fixes the wire spec, which is the
// actual defect.
class ShortedTurnsException : public CoilException {
public:
    explicit ShortedTurnsException(const std::string& message)
        : CoilException(ErrorCode::COIL_SHORTED_TURNS, message) {}
};

// ABT #1290 (Alf, 2026-09-20): a corner tighter than the wire's own minimum bend radius is not a
// numerical problem and not a tolerance: the layout asks the wire to do something the wire is not
// qualified to do, so the part it describes does not exist. Raised only under real winding, where
// the drawn corner is a physical claim; the ideal 2D layout is a separate contract and never
// raises it. Named after the fault rather than the symptom, for the same reason ShortedTurns is:
// a reader told "bend below the standard's mandrel" goes and looks at the former's corner radius
// or at the wire, which is where the defect is.
class UnwindableBendException : public CoilException {
public:
    explicit UnwindableBendException(const std::string& message)
        : CoilException(ErrorCode::COIL_WINDING_ERROR, message) {}
};

// ============================================================================
// Material Exceptions
// ============================================================================

class MaterialException : public OpenMagneticsException {
public:
    MaterialException(ErrorCode code, const std::string& message, const std::string& context = "")
        : OpenMagneticsException(code, message, context) {}
};

class MaterialDataMissingException : public MaterialException {
public:
    explicit MaterialDataMissingException(const std::string& materialName, const std::string& missingData = "")
        : MaterialException(ErrorCode::MATERIAL_DATA_MISSING, 
                           "Material data missing for: " + materialName + 
                           (missingData.empty() ? "" : " (missing: " + missingData + ")"), 
                           materialName) {}
};

/**
 * @brief A loss model was asked for a frequency outside the span its coefficients were fitted on.
 *
 * The Steinmetz ranges in MAS are fits to measured points over [minimumFrequency,
 * maximumFrequency]. Evaluating one outside that span is an extrapolation of a power law the data
 * never constrained: a 1-5 MHz fit evaluated at 100 kHz gave 2-18 % of 3C95's loss for nine MHz
 * grades (ABT #1456). Callers that sweep frequency must use get_fitted_span first; advisers drop
 * such materials through MagneticFilterLossModelFrequencySpan instead of catching this.
 */
class MaterialFrequencyOutOfSpanException : public MaterialException {
    std::string _materialName;
    double _frequency;
    double _spanMinimum;
    double _spanMaximum;
public:
    MaterialFrequencyOutOfSpanException(const std::string& materialName, const std::string& model,
                                        double frequency, double spanMinimum, double spanMaximum)
        : MaterialException(ErrorCode::MATERIAL_FREQUENCY_OUT_OF_SPAN,
                            "Material " + materialName + ": " + model + " loss coefficients are fitted over " +
                            std::to_string(spanMinimum) + " Hz to " + std::to_string(spanMaximum) +
                            " Hz; " + std::to_string(frequency) + " Hz is " +
                            (frequency < spanMinimum ? "below" : "above") +
                            " that span and would be an extrapolation of the fit",
                            materialName),
          _materialName(materialName), _frequency(frequency), _spanMinimum(spanMinimum), _spanMaximum(spanMaximum) {}

    const std::string& material_name() const noexcept { return _materialName; }
    double frequency() const noexcept { return _frequency; }
    double span_minimum() const noexcept { return _spanMinimum; }
    double span_maximum() const noexcept { return _spanMaximum; }
};

/**
 * @brief Complex permeability was asked for at a frequency outside the material's tabulated range.
 *
 * The complex permeability of a material is a table of measured (or, for materials that only carry
 * a frequency-dependent initial permeability, derived) points. Outside that table there is no data:
 * holding the last point made every higher frequency read the same mu' and mu'' (stale K081 data
 * tabulated 1 Hz to 1 MHz returned its 1 MHz point at 350 MHz and put a common-mode choke's
 * resonance there). Nothing clamps.
 *
 * ComplexPermeability::get_tabulated_complex_permeability throws this outside the tabulated range
 * (the advisers and filters, which stay on the data). ComplexPermeability::get_complex_permeability
 * estimates the missing range instead (ComplexPermeabilityEstimator, ABT #1861) and throws only the
 * subclass ComplexPermeabilityEstimateUnavailableException, for the materials it cannot estimate.
 */
class ComplexPermeabilityFrequencyOutOfRangeException : public MaterialException {
    std::string _materialName;
    double _frequency;
    double _rangeMinimum;
    double _rangeMaximum;
public:
    ComplexPermeabilityFrequencyOutOfRangeException(const std::string& materialName, double frequency,
                                                    double rangeMinimum, double rangeMaximum,
                                                    const std::string& reason = "")
        : MaterialException(ErrorCode::MATERIAL_FREQUENCY_OUT_OF_SPAN,
                            "Material " + materialName + ": complex permeability is tabulated from " +
                            std::to_string(rangeMinimum) + " Hz to " + std::to_string(rangeMaximum) +
                            " Hz; " + std::to_string(frequency) + " Hz is " +
                            (frequency < rangeMinimum ? "below" : "above") +
                            " that range and there is no data there" + (reason.empty() ? "" : "; " + reason),
                            materialName),
          _materialName(materialName), _frequency(frequency), _rangeMinimum(rangeMinimum), _rangeMaximum(rangeMaximum) {}

    const std::string& material_name() const noexcept { return _materialName; }
    double frequency() const noexcept { return _frequency; }
    double range_minimum() const noexcept { return _rangeMinimum; }
    double range_maximum() const noexcept { return _rangeMaximum; }
};

/**
 * @brief Complex permeability outside the tabulated range, for a material the estimator cannot cover.
 *
 * ComplexPermeability::get_complex_permeability estimates a material's missing frequency range from
 * the other materials of its class (ComplexPermeabilityEstimator, ABT #1861). It cannot for: a class
 * other than FERRITE MN_ZN / NI_ZN (powders, MgZn, a proprietary or missing composition), a material
 * whose complex permeability is derived from its initial permeability (no tabulated MAS data), a
 * class with fewer than 3 donors, and a material whose own data next to the gap cannot be matched.
 * reason() says which; the message names the material, the frequency and the tabulated range.
 */
class ComplexPermeabilityEstimateUnavailableException : public ComplexPermeabilityFrequencyOutOfRangeException {
    std::string _reason;
public:
    ComplexPermeabilityEstimateUnavailableException(const std::string& materialName, double frequency,
                                                    double rangeMinimum, double rangeMaximum, const std::string& reason)
        : ComplexPermeabilityFrequencyOutOfRangeException(materialName, frequency, rangeMinimum, rangeMaximum,
                                                          "and it cannot be estimated: " + reason),
          _reason(reason) {}

    const std::string& reason() const noexcept { return _reason; }
};

/**
 * @brief None of a requirement's frequencies lies inside the material's tabulated mu(f) range.
 *
 * A minimum-impedance requirement is judged only at the frequencies where the core material has
 * complex-permeability data; the points outside are not judged at all. When that leaves no point,
 * nothing about the requirement can be said for the part: it is not ranked, and the adviser reports
 * it with this error (MagneticAdviser::get_failed_candidates).
 */
class RequirementOutsideMaterialDataException : public MaterialException {
public:
    RequirementOutsideMaterialDataException(const std::string& materialName, const std::string& requirement,
                                            double requirementMinimum, double requirementMaximum,
                                            double rangeMinimum, double rangeMaximum,
                                            std::optional<std::pair<double, double>> measuredRange = std::nullopt)
        : MaterialException(ErrorCode::MATERIAL_FREQUENCY_OUT_OF_SPAN,
                            "Material " + materialName + ": complex permeability is tabulated from " +
                            std::to_string(rangeMinimum) + " Hz to " + std::to_string(rangeMaximum) +
                            " Hz" +
                            (measuredRange ? "; the part's measured impedance spans " + std::to_string(measuredRange->first) +
                                             " Hz to " + std::to_string(measuredRange->second) + " Hz"
                                           : std::string()) +
                            "; no " + requirement + " frequency (" + std::to_string(requirementMinimum) +
                            " Hz to " + std::to_string(requirementMaximum) +
                            " Hz) lies inside " + (measuredRange ? "either range" : "that range") +
                            ", so the requirement cannot be judged on this part",
                            materialName) {}
};

// A ferrite above its Curie point is paramagnetic: no loss model describes it, and the magnetic
// no longer works as designed (ABT #1485).
class MaterialAboveCurieTemperatureException : public MaterialException {
    std::string _materialName;
    double _temperature;
    double _curieTemperature;
public:
    MaterialAboveCurieTemperatureException(const std::string& materialName, double temperature, double curieTemperature)
        : MaterialException(ErrorCode::MATERIAL_ABOVE_CURIE_TEMPERATURE,
                            "Material " + materialName + ": core losses asked at " + std::to_string(temperature) +
                            " C, at or above its Curie temperature of " + std::to_string(curieTemperature) +
                            " C, where the material is no longer ferromagnetic",
                            materialName),
          _materialName(materialName), _temperature(temperature), _curieTemperature(curieTemperature) {}

    const std::string& material_name() const noexcept { return _materialName; }
    double temperature() const noexcept { return _temperature; }
    double curie_temperature() const noexcept { return _curieTemperature; }
};

// The initial permeability of a material was asked at a temperature below its Curie point but outside
// the temperatures its mu_i(T) table covers. There is no data there: the spline through the table is
// not a model of the material outside it (it ran 3C96's 180 C table to the Curie point and floored
// the result at 1). Above the Curie point the material is paramagnetic and mu_i = 1 is returned
// instead; Settings::allowMaterialDataExtrapolation turns this into an evaluation with a WARNING
// (ABT #1739).
class MaterialTemperatureOutOfRangeException : public MaterialException {
    std::string _materialName;
    double _temperature;
    double _tableMinimum;
    double _tableMaximum;
public:
    MaterialTemperatureOutOfRangeException(const std::string& materialName, double temperature,
                                           double tableMinimum, double tableMaximum,
                                           std::optional<double> curieTemperature)
        : MaterialException(ErrorCode::MATERIAL_TEMPERATURE_OUT_OF_RANGE,
                            "Material " + materialName + ": initial permeability is tabulated from " +
                            std::to_string(tableMinimum) + " C to " + std::to_string(tableMaximum) + " C" +
                            (curieTemperature ? ", Curie temperature " + std::to_string(curieTemperature.value()) + " C"
                                              : std::string(", no Curie temperature")) +
                            "; " + std::to_string(temperature) + " C is " +
                            (temperature < tableMinimum ? "below" : "above") +
                            " that table and there is no data there",
                            materialName),
          _materialName(materialName), _temperature(temperature), _tableMinimum(tableMinimum), _tableMaximum(tableMaximum) {}

    const std::string& material_name() const noexcept { return _materialName; }
    double temperature() const noexcept { return _temperature; }
    double table_minimum() const noexcept { return _tableMinimum; }
    double table_maximum() const noexcept { return _tableMaximum; }
};

// A caller explicitly asked for one core-loss model and the material cannot be evaluated with it
// (no data for it). Only an explicit request throws this; the default preference walks its cascade
// to the material's own model instead (ABT #1497).
// A loss factor tan(delta)/mu_i is measured at small signal (B <= 0.25 mT, IEC 60401-3). At a larger flux the
// hysteresis loss it leaves out dominates, so a loss computed from it there is not the material's loss (ABT #2006).
class LossFactorAboveSmallSignalFluxException : public MaterialException {
    std::string _materialName;
    double _magneticFluxDensityAmplitude;
    double _maximumMagneticFluxDensity;
public:
    LossFactorAboveSmallSignalFluxException(const std::string& materialName, double magneticFluxDensityAmplitude,
                                            double maximumMagneticFluxDensity)
        : MaterialException(ErrorCode::MATERIAL_LOSS_FACTOR_ABOVE_SMALL_SIGNAL_FLUX,
                            "Material " + materialName + ": its loss factor (tan(delta)/mu_i) is small-signal data, defined up to " +
                            std::to_string(maximumMagneticFluxDensity) + " T, and the core losses were asked at a flux amplitude of " +
                            std::to_string(magneticFluxDensityAmplitude) + " T, where it does not give the material's loss",
                            materialName),
          _materialName(materialName), _magneticFluxDensityAmplitude(magneticFluxDensityAmplitude),
          _maximumMagneticFluxDensity(maximumMagneticFluxDensity) {}

    const std::string& material_name() const noexcept { return _materialName; }
    double magnetic_flux_density_amplitude() const noexcept { return _magneticFluxDensityAmplitude; }
    double maximum_magnetic_flux_density() const noexcept { return _maximumMagneticFluxDensity; }
};

class RequestedCoreLossesModelNotAvailableException : public MaterialException {
    std::string _materialName;
    std::string _modelName;
    std::string _reason;
public:
    RequestedCoreLossesModelNotAvailableException(const std::string& materialName, const std::string& modelName,
                                                  const std::string& reason)
        : MaterialException(ErrorCode::MATERIAL_REQUESTED_LOSS_MODEL_NOT_AVAILABLE,
                            "Material " + materialName + ": core-loss model " + modelName +
                            " was explicitly requested but cannot be used: " + reason,
                            materialName),
          _materialName(materialName), _modelName(modelName), _reason(reason) {}

    const std::string& material_name() const noexcept { return _materialName; }
    const std::string& model_name() const noexcept { return _modelName; }
    const std::string& reason() const noexcept { return _reason; }
};

// ============================================================================
// Calculation Exceptions
// ============================================================================

class CalculationException : public OpenMagneticsException {
public:
    CalculationException(ErrorCode code, const std::string& message, const std::string& context = "")
        : OpenMagneticsException(code, message, context) {}
};

/**
 * @brief A material's magnetisation curve ends (at 1000 times its saturation field strength) below a flux
 * density that is still under its B_sat: the curve and the saturation data contradict each other there, and
 * no field strength can be read for that flux (InitialPermeability::get_magnetomotive_force_for_flux_density,
 * ABT #1942). Seen near the Curie point, where H_sat(T) goes to zero. The loss-temperature solve turns it, at a
 * heated iterate, into that operating point's NoThermalOperatingPointException (Alf, 2026-10-19).
 */
class MagnetisationCurveExhaustedException : public CalculationException {
public:
    explicit MagnetisationCurveExhaustedException(const std::string& message)
        : CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, message) {}
};

class NaNResultException : public CalculationException {
public:
    explicit NaNResultException(const std::string& calculationName)
        : CalculationException(ErrorCode::CALCULATION_NAN_RESULT, 
                              calculationName + " produced NaN result", calculationName) {}
};

class InvalidInputException : public CalculationException {
public:
    explicit InvalidInputException(const std::string& message)
        : CalculationException(ErrorCode::CALCULATION_INVALID_INPUT, message) {}
    
    InvalidInputException(ErrorCode code, const std::string& message)
        : CalculationException(code, message) {}
};

// The core saturates at the given operating point: a physics verdict on the design (this material,
// turn count, gap and core cannot carry that flux), not malformed input and not a failed calculation.
// An adviser rejects the candidate on exactly this exception; every other one still stops the run
// (ABT #1534, #1740). The message names the operating point and the flux against B_sat.
class SaturationException : public CalculationException {
    std::optional<int64_t> _numberTurns;
public:
    explicit SaturationException(const std::string& message)
        : CalculationException(ErrorCode::CALCULATION_CORE_SATURATES, message) {}
    // The turns-from-inductance solve refuses at a definite turn count. It carries that count so a caller
    // that can re-gap from it (the core adviser's seeding, ABT #1875) reads it here, never from the text.
    SaturationException(const std::string& message, int64_t numberTurns)
        : CalculationException(ErrorCode::CALCULATION_CORE_SATURATES, message), _numberTurns(numberTurns) {}

    // The turn count the refused solution had; empty when the refusal is not about one turn count.
    std::optional<int64_t> number_turns() const noexcept { return _numberTurns; }
};

// The thermal network's fixed point did not converge within its iteration limit (Temperature::
// calculateTemperatures with requireConvergence): the losses it was given have no steady temperature on
// that network. Its own type so the loss <-> temperature solve can tell it from every other calculation
// error (ABT #1876): inside that solve it means the operating point has no thermal operating point.
class ThermalNetworkNotConvergedException : public CalculationException {
public:
    explicit ThermalNetworkNotConvergedException(const std::string& message)
        : CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, message) {}
};

// A requirement point outside the material's mu(f) range is judged on the part's own measured
// common-mode |Z| (datasheet impedancePoints). This says why that measured curve cannot be read:
// a point with a non-positive or non-finite frequency or magnitude, two values at one frequency,
// or more than one zero-bias curve to choose from. Malformed data, never a reason to guess.
class InvalidDatasheetImpedanceException : public InvalidInputException {
public:
    InvalidDatasheetImpedanceException(const std::string& reference, const std::string& reason)
        : InvalidInputException(ErrorCode::INVALID_INPUT,
                                "Magnetic " + reference + ": datasheet commonModeChoke impedancePoints cannot be used: " + reason) {}
};

// ============================================================================
// Gap Exceptions
// ============================================================================

class GapException : public OpenMagneticsException {
public:
    GapException(ErrorCode code, const std::string& message, const std::string& context = "")
        : OpenMagneticsException(code, message, context) {}
    
    explicit GapException(const std::string& message)
        : OpenMagneticsException(ErrorCode::GAP_INVALID_DIMENSIONS, message, "") {}
};

class GapNotProcessedException : public GapException {
public:
    explicit GapNotProcessedException(const std::string& missingField)
        : GapException(ErrorCode::GAP_NOT_PROCESSED, 
                      "Gap " + missingField + " is not set") {}
};

// ============================================================================
// I/O Exceptions
// ============================================================================

class IOException : public OpenMagneticsException {
public:
    IOException(ErrorCode code, const std::string& message, const std::string& context = "")
        : OpenMagneticsException(code, message, context) {}
};

class FileNotFoundException : public IOException {
public:
    explicit FileNotFoundException(const std::string& filePath)
        : IOException(ErrorCode::IO_FILE_NOT_FOUND, 
                     "File not found: " + filePath, filePath) {}
};

class ParseException : public IOException {
public:
    explicit ParseException(const std::string& message, const std::string& filePath = "")
        : IOException(ErrorCode::IO_PARSE_ERROR, message, filePath) {}
};

// ============================================================================
// Model Exceptions  
// ============================================================================

class ModelException : public OpenMagneticsException {
public:
    ModelException(ErrorCode code, const std::string& message, const std::string& context = "")
        : OpenMagneticsException(code, message, context) {}
};

class ModelNotAvailableException : public ModelException {
public:
    explicit ModelNotAvailableException(const std::string& modelName, const std::string& materialName = "")
        : ModelException(ErrorCode::MODEL_NOT_AVAILABLE, 
                        "Model '" + modelName + "' not available" + 
                        (materialName.empty() ? "" : " for material: " + materialName)) {}
};

// ============================================================================
// Not Implemented Exception
// ============================================================================

class NotImplementedException : public OpenMagneticsException {
public:
    explicit NotImplementedException(const std::string& feature)
        : OpenMagneticsException(ErrorCode::NOT_IMPLEMENTED, 
                                feature + " is not implemented yet") {}
};

} // namespace OpenMagnetics
