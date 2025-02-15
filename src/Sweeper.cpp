#include "Sweeper.h"
#include "ComplexPermeability.h"
#include "Constants.h"
#include "WindingLosses.h"
#include "CoreLosses.h"
#include "Reluctance.h"
#include "InitialPermeability.h"
#include "Settings.h"
#include <cmath>


namespace OpenMagnetics {


Curve2D Sweeper::sweep_impedance_over_frequency(MagneticWrapper magnetic, double start, double stop, size_t numberElements, std::string mode, std::string title) {
    std::vector<double> frequencies;
    if (mode == "linear") {
        frequencies = linear_spaced_array(start, stop, numberElements);
    }
    else if (mode == "log") {
        frequencies = logarithmic_spaced_array(start, stop, numberElements);
    }
    else {
        throw std::runtime_error("Unknown spaced array mode");
    }

    std::vector<double> impedances;
    for (auto frequency : frequencies) {
        auto impedance = abs(OpenMagnetics::Impedance().calculate_impedance(magnetic, frequency));
        impedances.push_back(impedance);
    }

    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();
    return Curve2D(frequencies, impedances, title);
}

Curve2D Sweeper::sweep_winding_resistance_over_frequency(MagneticWrapper magnetic, double start, double stop, size_t numberElements, size_t windingIndex, double temperature, std::string mode, std::string title) {
    std::vector<double> frequencies;
    if (mode == "linear") {
        frequencies = linear_spaced_array(start, stop, numberElements);
    }
    else if (mode == "log") {
        frequencies = logarithmic_spaced_array(start, stop, numberElements);
    }
    else {
        throw std::runtime_error("Unknown spaced array mode");
    }

    auto magnetizingInductanceModel = MagnetizingInductance();
    auto turnsRatios = magnetic.get_mutable_coil().get_turns_ratios();
    auto magnetizingInductance = resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil()).get_magnetizing_inductance());
    double virtualCurrentRms = 1;
    std::vector<double> currentMask = {virtualCurrentRms * sqrt(2)};
    for (auto turnsRatio : turnsRatios) {
        currentMask.push_back(virtualCurrentRms * sqrt(2) * turnsRatio);
    }

    std::vector<double> effectiveResistances;
    for (auto frequency : frequencies) {
        auto operatingPoint = InputsWrapper::create_operating_point_with_sinusoidal_current_mask(frequency, magnetizingInductance, temperature, turnsRatios, currentMask);
        auto windingLossesPerWinding =  WindingLosses().calculate_losses(magnetic, operatingPoint, temperature).get_winding_losses_per_winding().value();

        auto proximityLossesPerharmonic = windingLossesPerWinding[windingIndex].get_proximity_effect_losses()->get_losses_per_harmonic();
        double proximityLosses = std::accumulate(proximityLossesPerharmonic.begin(), proximityLossesPerharmonic.end(), 0.0);
        auto skinLossesPerharmonic = windingLossesPerWinding[windingIndex].get_skin_effect_losses()->get_losses_per_harmonic();
        double skinLosses = std::accumulate(skinLossesPerharmonic.begin(), skinLossesPerharmonic.end(), 0.0);
        double lossesThisWinding = windingLossesPerWinding[windingIndex].get_ohmic_losses()->get_losses() + proximityLosses + skinLosses;

        double effectiveResistance = lossesThisWinding / pow(currentMask[windingIndex], 2);
        effectiveResistances.push_back(effectiveResistance);
        // std::cout << "frequency: " << frequency << std::endl;
        // std::cout << "effectiveResistance: " << effectiveResistance << std::endl;
    }

    return Curve2D(frequencies, effectiveResistances, title);
}

Curve2D Sweeper::sweep_resistance_over_frequency(MagneticWrapper magnetic, double start, double stop, size_t numberElements, double temperature, std::string mode, std::string title) {
    std::vector<double> frequencies;
    if (mode == "linear") {
        frequencies = linear_spaced_array(start, stop, numberElements);
    }
    else if (mode == "log") {
        frequencies = logarithmic_spaced_array(start, stop, numberElements);
    }
    else {
        throw std::runtime_error("Unknown spaced array mode");
    }

    auto magnetizingInductanceModel = MagnetizingInductance();
    auto turnsRatios = magnetic.get_mutable_coil().get_turns_ratios();
    auto magnetizingInductance = resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil()).get_magnetizing_inductance());
    double virtualCurrentRms = 1;
    std::vector<double> currentMask = {virtualCurrentRms * sqrt(2)};
    for (auto turnsRatio : turnsRatios) {
        currentMask.push_back(virtualCurrentRms * sqrt(2) * turnsRatio);
    }

    std::vector<double> effectiveResistances;
    for (auto frequency : frequencies) {
        auto operatingPoint = InputsWrapper::create_operating_point_with_sinusoidal_current_mask(frequency, magnetizingInductance, temperature, turnsRatios, currentMask);
        auto windingLosses =  WindingLosses().calculate_losses(magnetic, operatingPoint, temperature).get_winding_losses();

        double effectiveResistance = windingLosses / pow(currentMask[0], 2);
        effectiveResistances.push_back(effectiveResistance);
        // std::cout << "frequency: " << frequency << std::endl;
        // std::cout << "effectiveResistance: " << effectiveResistance << std::endl;
    }

    return Curve2D(frequencies, effectiveResistances, title);
}

Curve2D Sweeper::sweep_core_resistance_over_frequency(MagneticWrapper magnetic, double start, double stop, size_t numberElements, double temperature, std::string mode, std::string title) {
    std::vector<double> frequencies;
    if (mode == "linear") {
        frequencies = linear_spaced_array(start, stop, numberElements);
    }
    else if (mode == "log") {
        frequencies = logarithmic_spaced_array(start, stop, numberElements);
    }
    else {
        throw std::runtime_error("Unknown spaced array mode");
    }
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();

    auto magnetizingInductanceModel = MagnetizingInductance();
    auto magnetizingInductance = resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(core, coil).get_magnetizing_inductance());

    std::vector<double> coreResistances;
    auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(CoreLossesModels::STEINMETZ);
    for (auto frequency : frequencies) {
        auto coreResistance =  coreLossesModel->get_core_losses_series_resistance(core, frequency, temperature, magnetizingInductance);
        coreResistances.push_back(coreResistance);
    }

    return Curve2D(frequencies, coreResistances, title);
}

Curve2D Sweeper::sweep_core_losses_over_frequency(MagneticWrapper magnetic, OperatingPoint operatingPoint, double start, double stop, size_t numberElements, double temperature, std::string mode, std::string title) {
    std::vector<double> frequencies;
    if (mode == "linear") {
        frequencies = linear_spaced_array(start, stop, numberElements);
    }
    else if (mode == "log") {
        frequencies = logarithmic_spaced_array(start, stop, numberElements);
    }
    else {
        throw std::runtime_error("Unknown spaced array mode");
    }
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();

    auto magnetizingInductanceModel = MagnetizingInductance();
    auto magnetizingInductance = resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(core, coil).get_magnetizing_inductance());

    double numberWindings = coil.get_functional_description().size();
    double numberTurnsPrimary = coil.get_functional_description()[0].get_number_turns();
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();

    OpenMagnetics::InitialPermeability initialPermeabilityModel;
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory();
    auto initialPermeability = initialPermeabilityModel.get_initial_permeability(core.resolve_material(), temperature);
    auto magnetizingInductanceOutput = reluctanceModel->get_core_reluctance(core, initialPermeability);
    auto totalReluctance = magnetizingInductanceOutput.get_core_reluctance();


    std::vector<double> coreLossesPerFrequency;
    auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(CoreLossesModels::STEINMETZ);
    for (auto frequency : frequencies) {
        InputsWrapper::scale_time_to_frequency(operatingPoint, frequency);
        // operatingPoint = InputsWrapper::process_operating_point(operatingPoint, magnetizingInductance);
        OperatingPointExcitation excitation = InputsWrapper::get_primary_excitation(operatingPoint);

        if (numberWindings == 1 && excitation.get_current()) {
            InputsWrapper::set_current_as_magnetizing_current(&operatingPoint);
        }
        else if (InputsWrapper::is_multiport_inductor(operatingPoint)) {
            auto magnetizingCurrent = InputsWrapper::get_multiport_inductor_magnetizing_current(operatingPoint);
            excitation.set_magnetizing_current(magnetizingCurrent);
            operatingPoint.get_mutable_excitations_per_winding()[0] = excitation;
        }
        else if (excitation.get_voltage()) {
            auto voltage = operatingPoint.get_mutable_excitations_per_winding()[0].get_voltage().value();
            auto sampledVoltageWaveform = InputsWrapper::calculate_sampled_waveform(voltage.get_waveform().value(), frequency);

            auto magnetizingCurrent = InputsWrapper::calculate_magnetizing_current(excitation,
                                                                                   sampledVoltageWaveform,
                                                                                   magnetizingInductance,
                                                                                   false);

            auto sampledMagnetizingCurrentWaveform = InputsWrapper::calculate_sampled_waveform(magnetizingCurrent.get_waveform().value(), excitation.get_frequency());
            magnetizingCurrent.set_harmonics(InputsWrapper::calculate_harmonics_data(sampledMagnetizingCurrentWaveform, excitation.get_frequency()));
            magnetizingCurrent.set_processed(InputsWrapper::calculate_processed_data(magnetizingCurrent, sampledMagnetizingCurrentWaveform, false));

            excitation.set_magnetizing_current(magnetizingCurrent);
            operatingPoint.get_mutable_excitations_per_winding()[0] = excitation;
        }

        auto magneticFlux = OpenMagnetics::MagneticField::calculate_magnetic_flux(operatingPoint.get_mutable_excitations_per_winding()[0].get_magnetizing_current().value(), totalReluctance, numberTurnsPrimary);
        auto magneticFluxDensity = OpenMagnetics::MagneticField::calculate_magnetic_flux_density(magneticFlux, effectiveArea);
    
        excitation.set_magnetic_flux_density(magneticFluxDensity);

        auto coreLosses =  coreLossesModel->get_core_losses(core, excitation, temperature).get_core_losses();
        coreLossesPerFrequency.push_back(coreLosses);
    }

    return Curve2D(frequencies, coreLossesPerFrequency, title);
}

Curve2D Sweeper::sweep_winding_losses_over_frequency(MagneticWrapper magnetic, OperatingPoint operatingPoint, double start, double stop, size_t numberElements, double temperature, std::string mode, std::string title) {
    std::vector<double> frequencies;
    if (mode == "linear") {
        frequencies = linear_spaced_array(start, stop, numberElements);
    }
    else if (mode == "log") {
        frequencies = logarithmic_spaced_array(start, stop, numberElements);
    }
    else {
        throw std::runtime_error("Unknown spaced array mode");
    }

    auto magnetizingInductanceModel = MagnetizingInductance();
    auto turnsRatios = magnetic.get_mutable_coil().get_turns_ratios();
    auto magnetizingInductance = resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil()).get_magnetizing_inductance());
    double virtualCurrentRms = 1;
    std::vector<double> currentMask = {virtualCurrentRms * sqrt(2)};
    for (auto turnsRatio : turnsRatios) {
        currentMask.push_back(virtualCurrentRms * sqrt(2) * turnsRatio);
    }

    std::vector<double> windingLossesPerFrequency;
    for (auto frequency : frequencies) {
        InputsWrapper::scale_time_to_frequency(operatingPoint, frequency);
        operatingPoint = InputsWrapper::process_operating_point(operatingPoint, magnetizingInductance);

        auto windingLosses =  WindingLosses().calculate_losses(magnetic, operatingPoint, temperature).get_winding_losses();

        windingLossesPerFrequency.push_back(windingLosses);
        // std::cout << "frequency: " << frequency << std::endl;
        // std::cout << "effectiveResistance: " << effectiveResistance << std::endl;
    }

    return Curve2D(frequencies, windingLossesPerFrequency, title);
}


} // namespace OpenMagnetics
