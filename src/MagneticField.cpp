#include "MagneticField.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>

namespace OpenMagnetics {

ElectromagneticParameter MagneticField::get_magnetic_flux(ElectromagneticParameter magnetizingCurrent,
                                                          double reluctance,
                                                          double numberTurns,
                                                          double frequency) {
    ElectromagneticParameter magneticFlux;
    Waveform magneticFluxWaveform;
    std::vector<double> magneticFluxData;
    auto magnetizingCurrentWaveform = magnetizingCurrent.get_waveform().value();
    for (auto& datum : magnetizingCurrentWaveform.get_data()) {
        magneticFluxData.push_back(datum * numberTurns / reluctance);
    }

    if (magnetizingCurrentWaveform.get_time()) {
        magneticFluxWaveform.set_time(magnetizingCurrentWaveform.get_time());
    }

    magneticFluxWaveform.set_data(magneticFluxData);
    magneticFlux.set_waveform(magneticFluxWaveform);
    magneticFlux.set_harmonics(InputsWrapper::get_harmonics_data(magneticFluxWaveform, frequency));
    magneticFlux.set_processed(InputsWrapper::get_processed_data(magneticFlux, magneticFluxWaveform, true));

    return magneticFlux;
}
ElectromagneticParameter MagneticField::get_magnetic_flux_density(ElectromagneticParameter magneticFlux,
                                                                  double area,
                                                                  double frequency) {
    ElectromagneticParameter magneticFluxDensity;
    Waveform magneticFluxDensityWaveform;
    std::vector<double> magneticFluxDensityData;
    auto magneticFluxWaveform = magneticFlux.get_waveform().value();

    if (magneticFluxWaveform.get_time()) {
        magneticFluxDensityWaveform.set_time(magneticFluxWaveform.get_time());
    }

    for (auto& datum : magneticFluxWaveform.get_data()) {
        magneticFluxDensityData.push_back(datum / area);
    }

    magneticFluxDensityWaveform.set_data(magneticFluxDensityData);
    magneticFluxDensity.set_waveform(magneticFluxDensityWaveform);
    magneticFluxDensity.set_harmonics(InputsWrapper::get_harmonics_data(magneticFluxDensityWaveform, frequency));
    magneticFluxDensity.set_processed(
        InputsWrapper::get_processed_data(magneticFluxDensity, magneticFluxDensityWaveform, true));

    return magneticFluxDensity;
}

ElectromagneticParameter MagneticField::get_magnetic_field_strength(ElectromagneticParameter magneticFluxDensity,
                                                                    double initialPermeability,
                                                                    double frequency) {
    ElectromagneticParameter magneticFieldStrength;
    Waveform magneticFieldStrengthWaveform;
    std::vector<double> magneticFieldStrengthData;
    auto constants = Constants();
    auto magneticFluxDensityWaveform = magneticFluxDensity.get_waveform().value();

    if (magneticFluxDensityWaveform.get_time()) {
        magneticFieldStrengthWaveform.set_time(magneticFluxDensityWaveform.get_time());
    }

    for (auto& datum : magneticFluxDensityWaveform.get_data()) {
        magneticFieldStrengthData.push_back(datum / (initialPermeability * constants.vacuum_permeability));
    }

    magneticFieldStrengthWaveform.set_data(magneticFieldStrengthData);
    magneticFieldStrength.set_waveform(magneticFieldStrengthWaveform);
    magneticFieldStrength.set_harmonics(InputsWrapper::get_harmonics_data(magneticFieldStrengthWaveform, frequency));
    magneticFieldStrength.set_processed(
        InputsWrapper::get_processed_data(magneticFieldStrength, magneticFieldStrengthWaveform, true));

    return magneticFieldStrength;
}

} // namespace OpenMagnetics
