#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "support/Utils.h"
#include "constructive_models/Magnetic.h"

using namespace MAS;

namespace OpenMagnetics {


/**
 * @brief Structure holding topology-level waveforms for converter validation
 *
 * These waveforms represent the input/output signals of the converter
 * (input voltage, input current, output voltages, output currents),
 * NOT the magnetic component winding signals (which are OperatingPoints).
 */
class ConverterWaveforms {
public:
    ConverterWaveforms() = default;
    virtual ~ConverterWaveforms() = default;

private:
    Waveform inputVoltage;
    Waveform inputCurrent;
    std::vector<Waveform> outputVoltages;
    std::vector<Waveform> outputCurrents;
    std::string operatingPointName;
    double switchingFrequency = 0;

public:
    const Waveform& get_input_voltage() const { return inputVoltage; }
    Waveform& get_mutable_input_voltage() { return inputVoltage; }
    void set_input_voltage(const Waveform& value) { this->inputVoltage = value; }

    const Waveform& get_input_current() const { return inputCurrent; }
    Waveform& get_mutable_input_current() { return inputCurrent; }
    void set_input_current(const Waveform& value) { this->inputCurrent = value; }

    const std::vector<Waveform>& get_output_voltages() const { return outputVoltages; }
    std::vector<Waveform>& get_mutable_output_voltages() { return outputVoltages; }
    void set_output_voltages(const std::vector<Waveform>& value) { this->outputVoltages = value; }

    const std::vector<Waveform>& get_output_currents() const { return outputCurrents; }
    std::vector<Waveform>& get_mutable_output_currents() { return outputCurrents; }
    void set_output_currents(const std::vector<Waveform>& value) { this->outputCurrents = value; }

    const std::string& get_operating_point_name() const { return operatingPointName; }
    std::string& get_mutable_operating_point_name() { return operatingPointName; }
    void set_operating_point_name(const std::string& value) { this->operatingPointName = value; }

    const double& get_switching_frequency() const { return switchingFrequency; }
    double& get_mutable_switching_frequency() { return switchingFrequency; }
    void set_switching_frequency(const double& value) { this->switchingFrequency = value; }
};

class Topology {
public:
    bool _assertErrors = false;


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


    Topology(const json& j);
    Topology() {
    };

    static OperatingPointExcitation complete_excitation(Waveform currentWaveform, Waveform voltageWaveform, double switchingFrequency, std::string name);
    virtual DesignRequirements process_design_requirements() = 0;
    virtual std::vector<OperatingPoint> process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) = 0;
    virtual bool run_checks(bool assert = false);
    Inputs process();

};

} // namespace OpenMagnetics
