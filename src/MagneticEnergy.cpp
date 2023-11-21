#include "MagneticEnergy.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>

namespace OpenMagnetics {


double MagneticEnergy::get_ungapped_core_maximum_magnetic_energy(CoreWrapper core, OperatingPoint* operatingPoint){
    auto constants = Constants();
    double temperature = operatingPoint->get_conditions().get_ambient_temperature(); // TODO: Use a future calculated temperature
    double magneticFluxDensitySaturation = core.get_magnetic_flux_density_saturation(temperature);
    OpenMagnetics::InitialPermeability initialPermeability;
    double initialPermeabilityValue;
    if (operatingPoint != nullptr) {
        auto frequency = operatingPoint->get_excitations_per_winding()[0].get_frequency();
        initialPermeabilityValue = initialPermeability.get_initial_permeability(
            core.get_functional_description().get_material(), &temperature, nullptr, &frequency);
    }
    else {
        initialPermeabilityValue =
            initialPermeability.get_initial_permeability(core.get_functional_description().get_material());
    }

    double effective_volume = core.get_processed_description()->get_effective_parameters().get_effective_volume();

    double energyStoredInCore = 0.5 / (constants.vacuumPermeability * initialPermeabilityValue) * effective_volume *
                                  pow(magneticFluxDensitySaturation, 2);

    return energyStoredInCore;



}

double MagneticEnergy::get_gap_maximum_magnetic_energy(CoreGap gapInfo, double magneticFluxDensitySaturation, double* fringing_factor){
    auto constants = Constants();
    auto gap_length = gapInfo.get_length();
    auto gap_area = *(gapInfo.get_area());
    double fringing_factor_value;
    if (fringing_factor) {
        fringing_factor_value = *fringing_factor;
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

double MagneticEnergy::get_core_maximum_magnetic_energy(CoreWrapper core, OperatingPoint* operatingPoint){
    double totalEnergy = 0;
    double temperature = operatingPoint->get_conditions().get_ambient_temperature(); // TODO: Use a future calculated temperature
    double magneticFluxDensitySaturation = core.get_magnetic_flux_density_saturation(temperature);

    totalEnergy = get_ungapped_core_maximum_magnetic_energy(core, operatingPoint);
    for (auto& gap_info : core.get_functional_description().get_gapping()) {
        auto gapEnergy = get_gap_maximum_magnetic_energy(gap_info, magneticFluxDensitySaturation);
        totalEnergy += gapEnergy;
    }

    return totalEnergy;

}

DimensionWithTolerance MagneticEnergy::required_magnetic_energy(InputsWrapper inputs){
    DimensionWithTolerance desiredMagnetizingInductance = inputs.get_design_requirements().get_magnetizing_inductance();
    double magnetizingCurrentPeak = 0;
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
        magnetizingCurrentPeak = std::max(magnetizingCurrentPeak, InputsWrapper::get_primary_excitation(inputs.get_operating_point(operatingPointIndex)).get_magnetizing_current().value().get_processed().value().get_peak().value());
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
