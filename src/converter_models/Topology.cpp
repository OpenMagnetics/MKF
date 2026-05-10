#include "converter_models/PushPull.h"
#include "converter_models/Topology.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/Utils.h"
#include <cfloat>
#include <magic_enum.hpp>
#include "support/Exceptions.h"

namespace OpenMagnetics {

    // ----------------------------------------------------------------
    // SpiceSimulationConfig — central registry of per-topology defaults.
    //
    // Every topology that generates an ngspice netlist registers its
    // tuned values here. The struct's in-class defaults reflect a
    // generic "single-switch hard-switched DC-DC" baseline; topologies
    // diverge from those for documented physical reasons (snubber τ vs
    // Fs range, output cap vs power level, solver tolerances vs
    // switching events per period, etc.).
    //
    // Adding a new topology: insert a new entry below and override
    // `Topology::topology_kind()` in the concrete class.
    // ----------------------------------------------------------------
    const std::map<MAS::Topologies, SpiceSimulationConfig>&
    spice_simulation_defaults() {
        static const std::map<MAS::Topologies, SpiceSimulationConfig> defaults = [] {
            std::map<MAS::Topologies, SpiceSimulationConfig> m;

            // Boost — single-switch hard-switched, Fs 200 kHz – 2.2 MHz
            // (driven by reference-design EVMs TPS61089EVM-742,
            // TPS61178EVM-792, LM5122EVM-1PH). At 400 kHz a 1 µs snubber
            // τ would consume 40 % of the switching period and visibly
            // distort the inductor-current waveform, so the snubber is
            // sized for τ = 10 ns (negligible at all supported Fs).
            {
                SpiceSimulationConfig boost;
                boost.snubR = 100.0;       boost.snubC = 100e-12;
                boost.outputCapacitance = 100e-6;
                m[MAS::Topologies::BOOST_CONVERTER] = boost;
            }

            // Dab — dual active bridge, 8 switches in two H-bridges +
            // transformer, 4-corner phase-shift modulation. The original
            // hand-tuned values from Dab::generate_ngspice_circuit are
            // mirrored here verbatim so the registry is the single
            // source of truth. Snubber 1 kΩ/1 nF (τ = 1 µs) is fine at
            // the 50–250 kHz range Dab operates in. Solver tolerances
            // are looser than Boost's because the 4-position switching
            // produces 8 hard-switching events per period instead of 2.
            {
                SpiceSimulationConfig dab;
                dab.swModelVT = 2.5;       dab.swModelVH = 0.8;
                dab.swModelRON = 0.01;     dab.swModelROFF = 1e6;
                dab.snubR = 1e3;           dab.snubC = 1e-9;
                dab.diodeIS = 1e-12;       dab.diodeRS = 0.05;
                dab.outputCapacitance = 47e-6;
                dab.relTol = 0.01;         dab.absTol = 1e-7;
                dab.vnTol = 1e-4;
                dab.itl1 = 500;            dab.itl4 = 500;
                dab.method = "GEAR";       dab.trTol = 7.0;
                m[MAS::Topologies::DUAL_ACTIVE_BRIDGE_CONVERTER] = dab;
            }

            return m;
        }();
        return defaults;
    }

    SpiceSimulationConfig get_default_spice_config(MAS::Topologies t) {
        const auto& m = spice_simulation_defaults();
        auto it = m.find(t);
        if (it == m.end()) {
            throw std::invalid_argument(
                std::string("No SpiceSimulationConfig registered for topology ") +
                std::string(magic_enum::enum_name(t)) +
                ". Add an entry in Topology.cpp::spice_simulation_defaults().");
        }
        return it->second;
    }

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
        auto currentProcessed = Inputs::calculate_processed_data(currentWaveform, switchingFrequency, true);
        auto sampledCurrentWaveform = OpenMagnetics::Inputs::calculate_sampled_waveform(currentWaveform, switchingFrequency);
        auto currentHarmonics = OpenMagnetics::Inputs::calculate_harmonics_data(sampledCurrentWaveform, switchingFrequency);
        // Store the resampled (power-of-2) waveform so downstream consumers
        // (MagnetizingInductance, harmonics derivation, FFT-based pipelines)
        // get the standardized waveform contract MKF expects. The raw
        // analytical waveform is not necessarily a power-of-2 length (e.g.
        // DAB/PSFB/PSHB use 2*N+1 = 513 samples, AHB has variable size due to
        // discontinuity duplicates), and the size-check gates in
        // MagnetizingInductance.cpp would otherwise throw.
        current.set_waveform(sampledCurrentWaveform);
        current.set_processed(currentProcessed);
        current.set_harmonics(currentHarmonics);
        excitation.set_current(current);
        SignalDescriptor voltage;
        auto voltageProcessed = Inputs::calculate_processed_data(voltageWaveform, switchingFrequency, true);
        auto sampledVoltageWaveform = OpenMagnetics::Inputs::calculate_sampled_waveform(voltageWaveform, switchingFrequency);
        auto voltageHarmonics = OpenMagnetics::Inputs::calculate_harmonics_data(sampledVoltageWaveform, switchingFrequency);
        voltage.set_waveform(sampledVoltageWaveform);
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
