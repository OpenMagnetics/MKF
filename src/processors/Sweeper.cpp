#include "processors/Sweeper.h"
#include "physical_models/ComplexPermeability.h"
#include "Constants.h"
#include "physical_models/WindingLosses.h"
#include "physical_models/CoreLosses.h"
#include "physical_models/Reluctance.h"
#include "physical_models/InitialPermeability.h"
#include "support/Settings.h"
#include <cmath>


namespace OpenMagnetics {


Curve2D Sweeper::sweep_impedance_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, std::string mode, std::string title) {
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

Curve2D Sweeper::sweep_magnetizing_inductance_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, double temperature, std::string mode, std::string title) {
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
    auto staticMagnetizingInductance = resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil()).get_magnetizing_inductance());
    double virtualCurrentRms = 1;
    std::vector<double> currentMask = {virtualCurrentRms * sqrt(2)};
    for (auto turnsRatio : turnsRatios) {
        currentMask.push_back(virtualCurrentRms * sqrt(2) * turnsRatio);
    }

    std::vector<double> magnetizingInductances;
    for (auto frequency : frequencies) {
        auto operatingPoint = Inputs::create_operating_point_with_sinusoidal_current_mask(frequency, staticMagnetizingInductance, temperature, turnsRatios, currentMask);
        double magnetizingInductance = resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil(), &operatingPoint).get_magnetizing_inductance());
        magnetizingInductances.push_back(magnetizingInductance);
    }

    return Curve2D(frequencies, magnetizingInductances, title);
}

Curve2D Sweeper::sweep_magnetizing_inductance_over_temperature(Magnetic magnetic, double start, double stop, size_t numberElements, double frequency, std::string mode, std::string title) {
    std::vector<double> temperatures;
    if (mode == "linear") {
        temperatures = linear_spaced_array(start, stop, numberElements);
    }
    else if (mode == "log") {
        temperatures = logarithmic_spaced_array(start, stop, numberElements);
    }
    else {
        throw std::runtime_error("Unknown spaced array mode");
    }

    auto magnetizingInductanceModel = MagnetizingInductance();
    auto turnsRatios = magnetic.get_mutable_coil().get_turns_ratios();
    auto staticMagnetizingInductance = resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil()).get_magnetizing_inductance());
    double virtualCurrentRms = 1;
    std::vector<double> currentMask = {virtualCurrentRms * sqrt(2)};
    for (auto turnsRatio : turnsRatios) {
        currentMask.push_back(virtualCurrentRms * sqrt(2) * turnsRatio);
    }

    std::vector<double> magnetizingInductances;
    for (auto temperature : temperatures) {
        auto operatingPoint = Inputs::create_operating_point_with_sinusoidal_current_mask(frequency, staticMagnetizingInductance, temperature, turnsRatios, currentMask);
        double magnetizingInductance = resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil(), &operatingPoint).get_magnetizing_inductance());
        magnetizingInductances.push_back(magnetizingInductance);
    }

    return Curve2D(temperatures, magnetizingInductances, title);
}

Curve2D Sweeper::sweep_magnetizing_inductance_over_dc_bias(Magnetic magnetic, double start, double stop, size_t numberElements, double temperature, std::string mode, std::string title) {
    std::vector<double> currentOffsets;
    if (mode == "linear") {
        currentOffsets = linear_spaced_array(start, stop, numberElements);
    }
    else if (mode == "log") {
        currentOffsets = logarithmic_spaced_array(start, stop, numberElements);
    }
    else {
        throw std::runtime_error("Unknown spaced array mode");
    }

    auto magnetizingInductanceModel = MagnetizingInductance();

    Coil inductorCoil;
    inductorCoil.set_functional_description({magnetic.get_coil().get_functional_description()[0]});
    auto turnsRatios = inductorCoil.get_turns_ratios();
    auto staticMagnetizingInductance = resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), inductorCoil).get_magnetizing_inductance());
    double virtualCurrentRms = 1;
    std::vector<double> currentMask = {virtualCurrentRms * sqrt(2)};
    for (auto turnsRatio : turnsRatios) {
        currentMask.push_back(virtualCurrentRms * sqrt(2) * turnsRatio);
    }

    std::vector<double> magnetizingInductances;
    for (auto currentOffset : currentOffsets) {
        auto operatingPoint = Inputs::create_operating_point_with_sinusoidal_current_mask(defaults.measurementFrequency, staticMagnetizingInductance, temperature, turnsRatios, currentMask, currentOffset);
        double magnetizingInductance = resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), inductorCoil, &operatingPoint).get_magnetizing_inductance());
        magnetizingInductances.push_back(magnetizingInductance);
    }

    return Curve2D(currentOffsets, magnetizingInductances, title);
}

Curve2D Sweeper::sweep_winding_resistance_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, size_t windingIndex, double temperature, std::string mode, std::string title) {
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

    std::vector<double> effectiveResistances;
    for (auto frequency : frequencies) {
        auto effectiveResistance = WindingLosses::calculate_effective_resistance_of_winding(magnetic, windingIndex, frequency, temperature);
        effectiveResistances.push_back(effectiveResistance);
    }

    return Curve2D(frequencies, effectiveResistances, title);
}

Curve2D Sweeper::sweep_resistance_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, double temperature, std::string mode, std::string title) {
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
        auto operatingPoint = Inputs::create_operating_point_with_sinusoidal_current_mask(frequency, magnetizingInductance, temperature, turnsRatios, currentMask);
        auto windingLosses = WindingLosses().calculate_losses(magnetic, operatingPoint, temperature).get_winding_losses();

        double effectiveResistance = windingLosses / pow(operatingPoint.get_excitations_per_winding()[0].get_current()->get_processed()->get_rms().value(), 2);
        effectiveResistances.push_back(effectiveResistance);
    }

    return Curve2D(frequencies, effectiveResistances, title);
}

Curve2D Sweeper::sweep_core_resistance_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, double temperature, std::string mode, std::string title) {
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

Curve2D Sweeper::sweep_core_losses_over_frequency(Magnetic magnetic, OperatingPoint operatingPoint, double start, double stop, size_t numberElements, double temperature, std::string mode, std::string title) {
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
    CoreLosses coreLossesModel;
    coreLossesModel.set_core_losses_model_name(CoreLossesModels::STEINMETZ);

    for (auto frequency : frequencies) {
        Inputs::scale_time_to_frequency(operatingPoint, frequency, true);
        // operatingPoint = Inputs::process_operating_point(operatingPoint, magnetizingInductance);
        OperatingPointExcitation excitation = Inputs::get_primary_excitation(operatingPoint);

        if (numberWindings == 1 && excitation.get_current()) {
            Inputs::set_current_as_magnetizing_current(&operatingPoint);
        }
        else if (Inputs::is_multiport_inductor(operatingPoint, coil.get_isolation_sides())) {
            auto magnetizingCurrent = Inputs::get_multiport_inductor_magnetizing_current(operatingPoint);
            excitation.set_magnetizing_current(magnetizingCurrent);
            operatingPoint.get_mutable_excitations_per_winding()[0] = excitation;
        }
        else if (excitation.get_voltage()) {
            auto voltage = operatingPoint.get_mutable_excitations_per_winding()[0].get_voltage().value();
            auto sampledVoltageWaveform = Inputs::calculate_sampled_waveform(voltage.get_waveform().value(), frequency);

            auto magnetizingCurrent = Inputs::calculate_magnetizing_current(excitation,
                                                                                   sampledVoltageWaveform,
                                                                                   magnetizingInductance,
                                                                                   false);

            auto sampledMagnetizingCurrentWaveform = Inputs::calculate_sampled_waveform(magnetizingCurrent.get_waveform().value(), excitation.get_frequency());
            magnetizingCurrent.set_harmonics(Inputs::calculate_harmonics_data(sampledMagnetizingCurrentWaveform, excitation.get_frequency()));
            magnetizingCurrent.set_processed(Inputs::calculate_processed_data(magnetizingCurrent, sampledMagnetizingCurrentWaveform, false));

            excitation.set_magnetizing_current(magnetizingCurrent);
            operatingPoint.get_mutable_excitations_per_winding()[0] = excitation;
        }

        auto magneticFlux = OpenMagnetics::MagneticField::calculate_magnetic_flux(operatingPoint.get_mutable_excitations_per_winding()[0].get_magnetizing_current().value(), totalReluctance, numberTurnsPrimary);
        auto magneticFluxDensity = OpenMagnetics::MagneticField::calculate_magnetic_flux_density(magneticFlux, effectiveArea);
    
        excitation.set_magnetic_flux_density(magneticFluxDensity);

        auto coreLosses =  coreLossesModel.calculate_core_losses(core, excitation, temperature).get_core_losses();
        coreLossesPerFrequency.push_back(coreLosses);
    }

    return Curve2D(frequencies, coreLossesPerFrequency, title);
}

Curve2D Sweeper::sweep_winding_losses_over_frequency(Magnetic magnetic, OperatingPoint operatingPoint, double start, double stop, size_t numberElements, double temperature, std::string mode, std::string title) {
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

    std::vector<double> windingLossesPerFrequency;
    for (auto frequency : frequencies) {
        Inputs::scale_time_to_frequency(operatingPoint, frequency, true);
        operatingPoint = Inputs::process_operating_point(operatingPoint, magnetizingInductance);

        auto windingLosses =  WindingLosses().calculate_losses(magnetic, operatingPoint, temperature).get_winding_losses();
        auto windingLossesPerWinding =  WindingLosses().calculate_losses(magnetic, operatingPoint, temperature).get_winding_losses_per_winding().value();

        auto proximityLossesPerharmonic = windingLossesPerWinding[0].get_proximity_effect_losses()->get_losses_per_harmonic();

        windingLossesPerFrequency.push_back(windingLosses);
    }

    return Curve2D(frequencies, windingLossesPerFrequency, title);
}


} // namespace OpenMagnetics
