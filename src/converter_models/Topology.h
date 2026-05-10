#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include <CAS.hpp>
#include "processors/Inputs.h"
#include "support/Utils.h"
#include "constructive_models/Magnetic.h"
#include <map>
#include <optional>
#include <stdexcept>
#include <string>

using namespace MAS;

namespace OpenMagnetics {

/**
 * @brief Migrate old-format operating-point JSON (singular outputVoltage/
 *        outputCurrent) to the new plural format (outputVoltages array /
 *        outputCurrents array). Call this on the raw JSON object BEFORE
 *        passing it to quicktype's from_json, so old JSON inputs that
 *        predate the B3 schema change still parse correctly.
 *
 * Rules:
 *   - If "outputVoltage" (number) exists and "outputVoltages" doesn't →
 *     create "outputVoltages" = [outputVoltage], remove "outputVoltage".
 *   - If "outputCurrent" (number) exists and "outputCurrents" doesn't →
 *     create "outputCurrents" = [outputCurrent], remove "outputCurrent".
 *   - If both exist, the plural takes precedence (no migration needed).
 *
 * This handles the buck/boost/forward migration; flyback, LLC, DAB, etc.
 * already used the plural form and are unaffected.
 */
inline void migrate_operating_point_json(json& j) {
    if (j.contains("operatingPoints") && j["operatingPoints"].is_array()) {
        for (auto& op : j["operatingPoints"]) {
            if (op.contains("outputVoltage") && !op.contains("outputVoltages")) {
                op["outputVoltages"] = json::array({op["outputVoltage"]});
                op.erase("outputVoltage");
            }
            if (op.contains("outputCurrent") && !op.contains("outputCurrents")) {
                op["outputCurrents"] = json::array({op["outputCurrent"]});
                op.erase("outputCurrent");
            }
        }
    }
}


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

enum class ExtraComponentsMode {
    IDEAL,
    REAL
};

/**
 * @brief Configuration knobs for ngspice circuit generation.
 *
 * Centralizes every value that previously lived as a magic number inside
 * each topology's `generate_ngspice_circuit()` (PWM levels, switch model
 * thresholds, snubber RC, diode model, output filter, solver options).
 *
 * Two override points:
 *   1. **Per-topology defaults** — registered in
 *      `Topology.cpp::spice_simulation_defaults()`. Every topology MUST
 *      have an entry; missing entries throw rather than silently fall
 *      back to the universal struct defaults below (per the
 *      "no-fallbacks" rule).
 *   2. **Per-instance overrides** — call `set_spice_config()` on a
 *      topology instance to override any subset of fields, e.g. for
 *      tests that want a smaller `outputCapacitance` to shorten the
 *      RC settling tail.
 *
 * Usage inside a topology's circuit generator:
 *   auto cfg = spice_config();
 *   circuit << ".model SW1 SW VT=" << cfg.swModelVT << ...;
 */
struct SpiceSimulationConfig {
    // ---- Drive / switch ----
    double pwmHigh   = 5.0;       // PULSE high-level [V]
    double pwmRise   = 10e-9;     // PULSE rise time [s]
    double pwmFall   = 10e-9;     // PULSE fall time [s]
    double swModelVT = 2.5;       // ngspice SW model threshold [V]
    double swModelVH = 0.5;       // ngspice SW model hysteresis [V]

    // ---- Snubber RC across each switch ----
    double snubR = 100.0;         // [Ω]
    double snubC = 100e-12;       // [F]

    // ---- Output diode model (.model DIDEAL D(...)) ----
    double diodeIS = 1e-14;       // saturation current [A]
    double diodeRS = 1e-6;        // series resistance [Ω]

    // ---- Output filter ----
    double outputCapacitance              = 100e-6;  // [F]
    double outputCapInitialChargeFraction = 1.0;     // IC = Vout * fraction

    // ---- Solver / transient ----
    int    samplesPerPeriod = 200;
    double relTol  = 1e-3;
    double absTol  = 1e-9;
    double vnTol   = 1e-6;
    int    itl1    = 1000;
    int    itl4    = 1000;
    std::string method = "GEAR";
    double trTol   = 7.0;
};

/**
 * @brief Central registry of per-topology default `SpiceSimulationConfig`s.
 *
 * Every concrete topology that emits ngspice netlists must have an entry
 * here. The registry is keyed by `MAS::Topologies` so the compiler
 * catches typos.
 */
const std::map<MAS::Topologies, SpiceSimulationConfig>& spice_simulation_defaults();

/**
 * @brief Look up the default SpiceSimulationConfig for a given topology.
 *        Throws if the topology is not registered (no silent fallback).
 */
SpiceSimulationConfig get_default_spice_config(MAS::Topologies t);

class Topology {
public:
    bool _assertErrors = false;

private:
    std::optional<SpiceSimulationConfig> _spiceOverride;

public:
    /**
     * @brief Identify which MAS topology enum this concrete subclass
     *        represents. Used by `spice_config()` to look up defaults.
     *        Not pure-virtual to avoid breaking the few topologies that
     *        currently don't generate SPICE; those will throw if their
     *        `spice_config()` is called.
     */
    virtual MAS::Topologies topology_kind() const {
        throw std::logic_error(
            "topology_kind() not overridden — this topology does not "
            "support ngspice circuit generation");
    }

    /**
     * @brief Get the active `SpiceSimulationConfig`: per-instance
     *        override if set, otherwise the registered default.
     */
    SpiceSimulationConfig spice_config() const {
        if (_spiceOverride.has_value()) return _spiceOverride.value();
        return get_default_spice_config(topology_kind());
    }

    void set_spice_config(SpiceSimulationConfig cfg) { _spiceOverride = std::move(cfg); }
    void clear_spice_config() { _spiceOverride.reset(); }
    bool has_spice_config_override() const { return _spiceOverride.has_value(); }


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

    virtual std::vector<std::variant<Inputs, CAS::Inputs>> get_extra_components_inputs(
        ExtraComponentsMode mode,
        std::optional<Magnetic> magnetic = std::nullopt);

};

inline void to_json(json& j, const ConverterWaveforms& cw) {
    j["switchingFrequency"] = cw.get_switching_frequency();
    j["operatingPointName"] = cw.get_operating_point_name();
    
    auto inputVoltageWf = cw.get_input_voltage();
    if (!inputVoltageWf.get_data().empty()) {
        j["inputVoltage"] = json();
        j["inputVoltage"]["data"] = inputVoltageWf.get_data();
        if (inputVoltageWf.get_time()) {
            j["inputVoltage"]["time"] = inputVoltageWf.get_time().value();
        }
    }
    
    auto inputCurrentWf = cw.get_input_current();
    if (!inputCurrentWf.get_data().empty()) {
        j["inputCurrent"] = json();
        j["inputCurrent"]["data"] = inputCurrentWf.get_data();
        if (inputCurrentWf.get_time()) {
            j["inputCurrent"]["time"] = inputCurrentWf.get_time().value();
        }
    }
    
    auto outputVoltages = cw.get_output_voltages();
    j["outputVoltages"] = json::array();
    for (const auto& wf : outputVoltages) {
        json wfJson;
        wfJson["data"] = wf.get_data();
        if (wf.get_time()) {
            wfJson["time"] = wf.get_time().value();
        }
        j["outputVoltages"].push_back(wfJson);
    }
    
    auto outputCurrents = cw.get_output_currents();
    j["outputCurrents"] = json::array();
    for (const auto& wf : outputCurrents) {
        json wfJson;
        wfJson["data"] = wf.get_data();
        if (wf.get_time()) {
            wfJson["time"] = wf.get_time().value();
        }
        j["outputCurrents"].push_back(wfJson);
    }
}

} // namespace OpenMagnetics
