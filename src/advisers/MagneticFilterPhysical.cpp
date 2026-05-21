#include "advisers/MagneticFilter.h"
#include "advisers/MagneticFilterInternal.h"
#include "physical_models/MagnetizingInductance.h"
#include "processors/Inputs.h"
#include "support/Exceptions.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace OpenMagnetics {

// Phase 5: extracted from MagneticFilter.cpp.
// This translation unit owns the saturation / current-density filters:
//   - MagneticFilterSaturation (transformer-vs-inductor B detection)
//   - MagneticFilterDcCurrentDensity
//   - MagneticFilterEffectiveCurrentDensity
// The is_energy_storing_topology() helper lives in MagneticFilterInternal.h
// and is shared with the (still-monolithic) MagneticFilter.cpp.

std::pair<bool, double> MagneticFilterSaturation::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;

    // Transformer vs Inductor Detection
    // ==================================
    // Transformers: B is determined by applied VOLTAGE (Faraday's Law)
    //   - B_peak = V_peak / (N × Ae × ω)
    //   - No permeability iteration needed since magnetizing current is ideally small
    //
    // Inductors/Energy-storing: B is determined by CURRENT and permeability
    //   - Uses iterative calculation accounting for permeability rolloff with DC bias
    //
    // Detection priority:
    // 1. Use topology if specified (most reliable)
    // 2. Fall back to inductance field heuristic (minimum-only = transformer)
    //
    auto topology = inputs->get_design_requirements().get_topology();
    bool isTransformer;
    if (topology.has_value()) {
        // Use topology-based detection
        isTransformer = !is_energy_storing_topology(topology);
    } else {
        // Legacy heuristic: minimum-only inductance = transformer
        isTransformer = inputs->get_design_requirements().get_magnetizing_inductance().get_minimum() &&
                         !inputs->get_design_requirements().get_magnetizing_inductance().get_nominal() &&
                         !inputs->get_design_requirements().get_magnetizing_inductance().get_maximum();
    }

    for (auto operatingPoint : inputs->get_operating_points()) {
        double magneticFluxDensityPeak;
        auto magneticFluxDensitySaturation = magnetic->get_mutable_core().get_magnetic_flux_density_saturation(
            operatingPoint.get_conditions().get_ambient_temperature());
        
        if (isTransformer) {
            // For transformers, calculate B from voltage using MagnetizingInductance method
            // This avoids the permeability iteration issue with current-based calculation
            auto excitation = Inputs::get_primary_excitation(operatingPoint);
            double frequency = excitation.get_frequency();
            double voltagePeak = 0;
            double actualInductance = 0;
            
            bool hasVoltage = excitation.get_voltage() && excitation.get_voltage()->get_processed() && 
                excitation.get_voltage()->get_processed()->get_peak();
            bool hasCurrentWaveform = excitation.get_current() && excitation.get_current()->get_waveform();
            
            if (hasVoltage) {
                voltagePeak = excitation.get_voltage()->get_processed()->get_peak().value();
            } else if (hasCurrentWaveform) {
                // Derive voltage from V = L·di/dt using the actual gapped-core inductance
                OpenMagnetics::MagnetizingInductance magnetizingInductanceCalc;
                auto inductanceOutput = magnetizingInductanceCalc.calculate_inductance_from_number_turns_and_gapping(
                    *magnetic, &operatingPoint);

                if (inductanceOutput.get_magnetizing_inductance().get_nominal()) {
                    actualInductance = inductanceOutput.get_magnetizing_inductance().get_nominal().value();
                    if (actualInductance > 0) {
                        auto inducedVoltage = Inputs::calculate_induced_voltage(excitation, actualInductance);
                        if (inducedVoltage.get_processed() && inducedVoltage.get_processed()->get_peak()) {
                            voltagePeak = inducedVoltage.get_processed()->get_peak().value();
                        }
                    }
                }
            }
            
            double numberTurns = magnetic->get_coil().get_functional_description()[0].get_number_turns();
            OpenMagnetics::MagnetizingInductance magnetizingInductanceObj;

            // For non-sinusoidal waveforms (e.g. DAB square wave), the sinusoidal formula
            // V_peak/(N·Ae·ω) underestimates B_peak by ~36%. Integrate V·dt instead.
            double maxVoltSeconds = 0.0;
            if (hasVoltage && excitation.get_voltage()->get_waveform() &&
                excitation.get_voltage()->get_waveform()->get_time()) {
                auto voltageWaveform = excitation.get_voltage()->get_waveform().value();
                const auto& data = voltageWaveform.get_data();
                auto time = voltageWaveform.get_time().value();
                double integral = 0.0;
                for (size_t j = 0; j + 1 < std::min(data.size(), time.size()); ++j) {
                    integral += data[j] * (time[j + 1] - time[j]);
                    maxVoltSeconds = std::max(maxVoltSeconds, std::abs(integral));
                }
            }

            double effectiveArea = magnetic->get_mutable_core().get_processed_description()->get_effective_parameters().get_effective_area();
            if (maxVoltSeconds > 0 && numberTurns > 0 && effectiveArea > 0) {
                magneticFluxDensityPeak = maxVoltSeconds / (numberTurns * effectiveArea);
            } else {
                magneticFluxDensityPeak = magnetizingInductanceObj.calculate_flux_density_peak_from_voltage(
                    magnetic->get_mutable_core(), numberTurns, voltagePeak, frequency);
            }
        } else {
            // For inductors/energy-storing converters, use the standard calculation
            OpenMagnetics::MagnetizingInductance magnetizingInductanceObj;
            auto magneticFluxDensity = magnetizingInductanceObj.calculate_inductance_and_magnetic_flux_density(
                magnetic->get_core(), magnetic->get_coil(), &operatingPoint).second;
            magneticFluxDensityPeak = magneticFluxDensity.get_processed().value().get_peak().value();
        }

        // Saturation scoring: dimensionless ratio of operating peak flux density to
        // material saturation. Smaller is better (more headroom). The previous
        // implementation used fabs(Bsat - Bpeak), which under min-max+invert
        // *rewards* candidates that operate close to saturation - the opposite of
        // engineering intent. Ratio Bpeak/Bsat preserves the "more headroom = higher
        // normalized score" semantics correctly.
        double bRatio = (magneticFluxDensitySaturation > 0)
            ? (magneticFluxDensityPeak / magneticFluxDensitySaturation)
            : 1.0;
        scoring += bRatio;

        bool isSaturated = magneticFluxDensityPeak > magneticFluxDensitySaturation;

        if (isSaturated) {
            return {false, 0.0};
        }
    }

    scoring /= inputs->get_operating_points().size();

    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterDcCurrentDensity::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;

    if (inputs->get_operating_points().size() > 0 && magnetic->get_mutable_coil().get_functional_description().size() != inputs->get_operating_points()[0].get_excitations_per_winding().size()) {
        return {false, 0.0};        
    }

    for (const auto& operatingPoint : inputs->get_operating_points()) {
        for (size_t windingIndex = 0; windingIndex < magnetic->get_mutable_coil().get_functional_description().size(); ++windingIndex) {
            auto excitation = operatingPoint.get_excitations_per_winding()[windingIndex];
            if (!excitation.get_current()) { 
                throw InvalidInputException(ErrorCode::MISSING_DATA, "Current is missing in excitation");
            }
            auto current = excitation.get_current().value();
            auto wire = magnetic->get_mutable_coil().resolve_wire(windingIndex);
            auto dcCurrentDensity = wire.calculate_dc_current_density(current) / magnetic->get_mutable_coil().get_functional_description()[windingIndex].get_number_parallels();

            // DC current density scoring: ratio of actual to maximum allowed.
            // Smaller is better (more headroom against thermal limit). See note on
            // MagneticFilterSaturation for why fabs(max - actual) was wrong.
            scoring += dcCurrentDensity / defaults.maximumCurrentDensity;
            if (dcCurrentDensity > defaults.maximumCurrentDensity) {
                return {false, 0.0};
            }
        }
    }

    scoring /= inputs->get_operating_points().size();

    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterEffectiveCurrentDensity::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;

    if (inputs->get_operating_points().size() > 0 && magnetic->get_mutable_coil().get_functional_description().size() != inputs->get_operating_points()[0].get_excitations_per_winding().size()) {
        return {false, 0.0};        
    }

    for (const auto& operatingPoint : inputs->get_operating_points()) {
        for (size_t windingIndex = 0; windingIndex < magnetic->get_mutable_coil().get_functional_description().size(); ++windingIndex) {
            auto excitation = operatingPoint.get_excitations_per_winding()[windingIndex];
            if (!excitation.get_current()) {
                throw InvalidInputException(ErrorCode::MISSING_DATA, "Current is missing in excitation");
            }
            auto current = excitation.get_current().value();
            auto wire = magnetic->get_mutable_coil().resolve_wire(windingIndex);
            auto effectiveCurrentDensity = wire.calculate_effective_current_density(current, operatingPoint.get_conditions().get_ambient_temperature()) / magnetic->get_mutable_coil().get_functional_description()[windingIndex].get_number_parallels();

            // Effective (AC+DC, temperature-corrected) current density scoring:
            // ratio of actual to maximum allowed. Smaller is better. See note on
            // MagneticFilterSaturation for the rationale.
            scoring += effectiveCurrentDensity / defaults.maximumEffectiveCurrentDensity;
            if (effectiveCurrentDensity > defaults.maximumEffectiveCurrentDensity) {
                return {false, 0.0};
            }
        }
    }

    scoring /= inputs->get_operating_points().size();

    return {valid, scoring};
}

} // namespace OpenMagnetics
