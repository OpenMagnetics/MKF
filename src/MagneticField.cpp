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

SignalDescriptor MagneticField::calculate_magnetic_flux(SignalDescriptor magnetizingCurrent,
                                                          double reluctance,
                                                          double numberTurns,
                                                          double frequency) {
    SignalDescriptor magneticFlux;
    Waveform magneticFluxWaveform;
    std::vector<double> magneticFluxData;
    auto magnetizingCurrentWaveform = magnetizingCurrent.get_waveform().value();

    if (InputsWrapper::is_waveform_sampled(magnetizingCurrentWaveform)) {
        if (magnetizingCurrentWaveform.get_data().size() > 0 && ((magnetizingCurrentWaveform.get_data().size() & (magnetizingCurrentWaveform.get_data().size() - 1)) != 0)) {
            std::cout << "magnetizingCurrentWaveform.get_data().size():" << magnetizingCurrentWaveform.get_data().size() << std::endl;
            throw std::invalid_argument("magnetizingCurrentWaveform vector size is not a power of 2");
        }
    }
    else {
         magnetizingCurrentWaveform = InputsWrapper::calculate_sampled_waveform(magnetizingCurrentWaveform, frequency);
    }

    for (auto& datum : magnetizingCurrentWaveform.get_data()) {
        magneticFluxData.push_back(datum * numberTurns / reluctance);
    }

    if (magnetizingCurrentWaveform.get_time()) {
        magneticFluxWaveform.set_time(magnetizingCurrentWaveform.get_time());
    }

    magneticFluxWaveform.set_data(magneticFluxData);
    magneticFlux.set_waveform(magneticFluxWaveform);
    // magneticFlux.set_harmonics(InputsWrapper::calculate_harmonics_data(magneticFluxWaveform, frequency));
    // magneticFlux.set_processed(InputsWrapper::calculate_processed_data(magneticFlux, magneticFluxWaveform, true, false));

    return magneticFlux;
}
SignalDescriptor MagneticField::calculate_magnetic_flux_density(SignalDescriptor magneticFlux,
                                                                  double area,
                                                                  double frequency) {
    SignalDescriptor magneticFluxDensity;
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
    magneticFluxDensity.set_harmonics(InputsWrapper::calculate_harmonics_data(magneticFluxDensityWaveform, frequency));
    magneticFluxDensity.set_processed(
        InputsWrapper::calculate_processed_data(magneticFluxDensity, magneticFluxDensityWaveform, true, false));

    return magneticFluxDensity;
}

SignalDescriptor MagneticField::calculate_magnetic_field_strength(SignalDescriptor magneticFluxDensity,
                                                                    double initialPermeability,
                                                                    double frequency) {
    SignalDescriptor magneticFieldStrength;
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
    magneticFieldStrength.set_harmonics(InputsWrapper::calculate_harmonics_data(magneticFieldStrengthWaveform, frequency));
    magneticFieldStrength.set_processed(
        InputsWrapper::calculate_processed_data(magneticFieldStrength, magneticFieldStrengthWaveform, true, false));

    return magneticFieldStrength;
}

} // namespace OpenMagnetics
