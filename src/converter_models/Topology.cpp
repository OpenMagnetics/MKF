#include "converter_models/PushPull.h"
#include "converter_models/Topology.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/Utils.h"
#include <cfloat>
#include "support/Exceptions.h"

namespace OpenMagnetics {

    Topology::Topology(const json& j) {
        // Base class constructor - derived classes should handle their own JSON deserialization
    }

    OperatingPointExcitation Topology::complete_excitation(Waveform currentWaveform, Waveform voltageWaveform, double switchingFrequency, std::string name) {
        if (switchingFrequency <= 0 || !std::isfinite(switchingFrequency)) {
            throw std::invalid_argument("complete_excitation: Invalid switchingFrequency: " + std::to_string(switchingFrequency));
        }

        OperatingPointExcitation excitation;
        excitation.set_frequency(switchingFrequency);
        
        SignalDescriptor current;
        current.set_waveform(currentWaveform);
        auto currentProcessed = Inputs::calculate_processed_data(currentWaveform, switchingFrequency, true);
        auto sampledCurrentWaveform = OpenMagnetics::Inputs::calculate_sampled_waveform(currentWaveform, switchingFrequency);
        auto currentHarmonics = OpenMagnetics::Inputs::calculate_harmonics_data(sampledCurrentWaveform, switchingFrequency);
        current.set_processed(currentProcessed);
        current.set_harmonics(currentHarmonics);
        excitation.set_current(current);
        SignalDescriptor voltage;
        voltage.set_waveform(voltageWaveform);
        auto voltageProcessed = Inputs::calculate_processed_data(voltageWaveform, switchingFrequency, true);
        auto sampledVoltageWaveform = OpenMagnetics::Inputs::calculate_sampled_waveform(voltageWaveform, switchingFrequency);
        auto voltageHarmonics = OpenMagnetics::Inputs::calculate_harmonics_data(sampledVoltageWaveform, switchingFrequency);
        voltage.set_processed(voltageProcessed);
        voltage.set_harmonics(voltageHarmonics);
        excitation.set_voltage(voltage);
        excitation.set_name(name);
        return excitation;
    }

    bool Topology::run_checks(bool assert) {
        throw NotImplementedException("Not implemented");
    }

    Inputs Topology::process() {
        run_checks(_assertErrors);

        Inputs inputs;
        auto designRequirements = process_design_requirements();
        std::vector<double> turnsRatios;
        for (auto turnsRatio : designRequirements.get_turns_ratios()) {
            turnsRatios.push_back(resolve_dimensional_values(turnsRatio));
        }
        auto desiredMagnetizingInductance = resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
        auto operatingPoints = process_operating_points(turnsRatios, desiredMagnetizingInductance);
        for (size_t operatingPointIndex = 0; operatingPointIndex < operatingPoints.size(); ++operatingPointIndex) {
            operatingPoints[operatingPointIndex] = Inputs::prune_harmonics(operatingPoints[operatingPointIndex], Defaults().harmonicAmplitudeThreshold);
        }

        inputs.set_design_requirements(designRequirements);
        inputs.set_operating_points(operatingPoints);

        return inputs;
    }


    std::vector<std::variant<Inputs, CAS::Inputs>> Topology::get_extra_components_inputs(
        ExtraComponentsMode mode,
        std::optional<Magnetic> magnetic)
    {
        if (mode == ExtraComponentsMode::REAL && !magnetic.has_value()) {
            throw std::invalid_argument("get_extra_components_inputs: mode REAL requires a designed magnetic");
        }
        return {};
    }

} // namespace OpenMagnetics
