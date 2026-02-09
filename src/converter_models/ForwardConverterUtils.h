#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"


namespace OpenMagnetics {
using namespace MAS;
namespace ForwardConverterUtils {

/**
 * @brief Collect input voltages from dimension with tolerance
 */
inline void collect_input_voltages(const DimensionWithTolerance& inputVoltage,
                                   std::vector<double>& voltages, 
                                   std::vector<std::string>& names) {
    voltages.clear();
    names.clear();
    
    if (inputVoltage.get_nominal()) {
        voltages.push_back(inputVoltage.get_nominal().value());
        names.push_back("Nom.");
    }
    if (inputVoltage.get_minimum()) {
        voltages.push_back(inputVoltage.get_minimum().value());
        names.push_back("Min.");
    }
    if (inputVoltage.get_maximum()) {
        voltages.push_back(inputVoltage.get_maximum().value());
        names.push_back("Max.");
    }
}

/**
 * @brief Create isolation sides vector for forward converters
 */
inline std::vector<IsolationSide> create_isolation_sides(size_t numSecondaries, bool hasDemagnetizationWinding) {
    std::vector<IsolationSide> isolationSides;
    
    // Primary
    isolationSides.push_back(get_isolation_side_from_index(0));
    
    // Optional demagnetization winding
    if (hasDemagnetizationWinding) {
        isolationSides.push_back(get_isolation_side_from_index(0));
    }
    
    // Secondaries
    for (size_t i = 0; i < numSecondaries; ++i) {
        isolationSides.push_back(get_isolation_side_from_index(i + 1));
    }
    
    return isolationSides;
}

/**
 * @brief Common run_checks implementation for forward converters
 */
template<typename Converter>
inline bool run_checks_common(Converter* converter, bool assert) {
    if (converter->get_operating_points().size() == 0) {
        if (!assert) {
            return false;
        }
        throw InvalidInputException(ErrorCode::MISSING_DATA, "At least one operating point is needed");
    }
    for (size_t i = 1; i < converter->get_operating_points().size(); ++i) {
        if (converter->get_operating_points()[i].get_output_voltages().size() != 
            converter->get_operating_points()[0].get_output_voltages().size()) {
            if (!assert) {
                return false;
            }
            throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, 
                "Different operating points cannot have different number of output voltages");
        }
        if (converter->get_operating_points()[i].get_output_currents().size() != 
            converter->get_operating_points()[0].get_output_currents().size()) {
            if (!assert) {
                return false;
            }
            throw InvalidInputException(ErrorCode::INVALID_DESIGN_REQUIREMENTS, 
                "Different operating points cannot have different number of output currents");
        }
    }
    if (!converter->get_input_voltage().get_nominal() && 
        !converter->get_input_voltage().get_maximum() && 
        !converter->get_input_voltage().get_minimum()) {
        if (!assert) {
            return false;
        }
        throw InvalidInputException(ErrorCode::MISSING_DATA, "No input voltage introduced");
    }

    return true;
}

} // namespace ForwardConverterUtils
} // namespace OpenMagnetics
