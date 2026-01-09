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
    
    // Material related errors (300-399)
    MATERIAL_NOT_FOUND = 300,
    MATERIAL_DATA_MISSING = 301,
    MATERIAL_INVALID_PROPERTY = 302,
    
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
    CALCULATION_TIMEOUT = 1301
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
        case ErrorCode::MATERIAL_NOT_FOUND: return "MATERIAL_NOT_FOUND";
        case ErrorCode::MATERIAL_DATA_MISSING: return "MATERIAL_DATA_MISSING";
        case ErrorCode::MATERIAL_INVALID_PROPERTY: return "MATERIAL_INVALID_PROPERTY";
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
        default: return "UNKNOWN";
    }
}

/**
 * @brief Base exception class for all OpenMagnetics exceptions
 */
class OpenMagneticsException : public std::exception {
protected:
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

// ============================================================================
// Calculation Exceptions
// ============================================================================

class CalculationException : public OpenMagneticsException {
public:
    CalculationException(ErrorCode code, const std::string& message, const std::string& context = "")
        : OpenMagneticsException(code, message, context) {}
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
