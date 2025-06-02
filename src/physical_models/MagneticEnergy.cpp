#include "physical_models/MagneticEnergy.h"
#include "Defaults.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>

namespace OpenMagnetics {


double MagneticEnergy::get_ungapped_core_maximum_magnetic_energy(Core core, std::optional<OperatingPoint> operatingPoint, bool saturationProportion){
    auto constants = Constants();
    double temperature = Defaults().ambientTemperature;
    if (operatingPoint) {
        temperature = operatingPoint->get_conditions().get_ambient_temperature(); // TODO: Use a future calculated temperature
    }
    if (operatingPoint) {
        auto frequency = operatingPoint->get_excitations_per_winding()[0].get_frequency();
        return get_ungapped_core_maximum_magnetic_energy(core, temperature, frequency, saturationProportion);
    }
    else {
        return get_ungapped_core_maximum_magnetic_energy(core, temperature, saturationProportion);
    }
}

double MagneticEnergy::get_ungapped_core_maximum_magnetic_energy(Core core, double temperature, std::optional<double> frequency, bool saturationProportion){
    auto constants = Constants();
    double magneticFluxDensitySaturation = core.get_magnetic_flux_density_saturation(temperature, saturationProportion);
    OpenMagnetics::InitialPermeability initialPermeability;
    auto coreMaterial = core.resolve_material();

    double initialPermeabilityValue = initialPermeability.get_initial_permeability(coreMaterial, temperature, std::nullopt, frequency);
    double effective_volume = core.get_processed_description()->get_effective_parameters().get_effective_volume();
    double energyStoredInCore = 0.5 / (constants.vacuumPermeability * initialPermeabilityValue) * effective_volume * pow(magneticFluxDensitySaturation, 2);

    return energyStoredInCore;
}

double MagneticEnergy::get_gap_maximum_magnetic_energy(CoreGap gapInfo, double magneticFluxDensitySaturation, std::optional<double> fringing_factor){
    auto constants = Constants();
    auto gap_length = gapInfo.get_length();
    auto gap_area = *(gapInfo.get_area());
    double fringing_factor_value;
    if (fringing_factor) {
        fringing_factor_value = fringing_factor.value();
    }
    else {
        auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(
            magic_enum::enum_cast<OpenMagnetics::ReluctanceModels>(_models["gapReluctance"]).value());

        fringing_factor_value = reluctanceModel->get_gap_reluctance(gapInfo).get_fringing_factor();
    }

    double energyStoredInGap = 0.5 / constants.vacuumPermeability * gap_length * gap_area * fringing_factor_value *
                                  pow(magneticFluxDensitySaturation, 2);

    return energyStoredInGap;

}

double MagneticEnergy::calculate_core_maximum_magnetic_energy(Core core, std::optional<OperatingPoint> operatingPoint, bool saturationProportion){
    double temperature = Defaults().ambientTemperature;
    if (operatingPoint) {
        temperature = operatingPoint->get_conditions().get_ambient_temperature(); // TODO: Use a future calculated temperature
    }
    return calculate_core_maximum_magnetic_energy(core, temperature, saturationProportion);
}

double MagneticEnergy::calculate_core_maximum_magnetic_energy(Core core, double temperature, std::optional<double> frequency, bool saturationProportion){
    double totalEnergy = 0;
    double magneticFluxDensitySaturation = core.get_magnetic_flux_density_saturation(temperature, saturationProportion);

    totalEnergy = get_ungapped_core_maximum_magnetic_energy(core, temperature, frequency);
    for (auto& gap_info : core.get_functional_description().get_gapping()) {
        auto gapEnergy = get_gap_maximum_magnetic_energy(gap_info, magneticFluxDensitySaturation);
        totalEnergy += gapEnergy;
    }

    return totalEnergy;
}

DimensionWithTolerance MagneticEnergy::calculate_required_magnetic_energy(Inputs inputs){
    DimensionWithTolerance desiredMagnetizingInductance = inputs.get_design_requirements().get_magnetizing_inductance();
    double magnetizingCurrentPeak = 0;
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
        if (!Inputs::get_primary_excitation(inputs.get_operating_point(operatingPointIndex)).get_magnetizing_current()) {
            throw std::runtime_error("Missing magnetizing current");
        }
        if (!Inputs::get_primary_excitation(inputs.get_operating_point(operatingPointIndex)).get_magnetizing_current()->get_processed()) {
            auto excitation = Inputs::get_primary_excitation(inputs.get_operating_point(operatingPointIndex));
            auto magnetizingCurrent = excitation.get_magnetizing_current().value();
            auto magnetizingCurrentExcitationWaveform = magnetizingCurrent.get_waveform().value();
            auto sampledCurrentWaveform = Inputs::calculate_sampled_waveform(magnetizingCurrentExcitationWaveform, excitation.get_frequency());
            magnetizingCurrent.set_harmonics(Inputs::calculate_harmonics_data(sampledCurrentWaveform, excitation.get_frequency()));
            magnetizingCurrent.set_processed(Inputs::calculate_processed_data(magnetizingCurrent, sampledCurrentWaveform, true, magnetizingCurrent.get_processed()));
            excitation.set_magnetizing_current(magnetizingCurrent);
            inputs.get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[0] = excitation;
        }
        magnetizingCurrentPeak = std::max(magnetizingCurrentPeak, Inputs::get_primary_excitation(inputs.get_operating_point(operatingPointIndex)).get_magnetizing_current().value().get_processed().value().get_peak().value());
    }
    DimensionWithTolerance magneticEnergyRequirement;
    auto get_energy = [magnetizingCurrentPeak](double magnetizingInductance)
    {
        return magnetizingInductance * pow(magnetizingCurrentPeak, 2) / 2;
    };

    if (desiredMagnetizingInductance.get_maximum()) {
        magneticEnergyRequirement.set_minimum(get_energy(*desiredMagnetizingInductance.get_maximum()));
    }
    if (desiredMagnetizingInductance.get_minimum()) {
        magneticEnergyRequirement.set_maximum(get_energy(*desiredMagnetizingInductance.get_minimum()));
    }
    if (desiredMagnetizingInductance.get_nominal()) {
        magneticEnergyRequirement.set_nominal(get_energy(*desiredMagnetizingInductance.get_nominal()));
    }

    return magneticEnergyRequirement;
}

} // namespace OpenMagnetics
