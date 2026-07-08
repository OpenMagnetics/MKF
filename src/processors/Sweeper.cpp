#include "processors/Sweeper.h"
#include "physical_models/ComplexPermeability.h"
#include "Constants.h"
#include "physical_models/WindingLosses.h"
#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/CoreLosses.h"
#include "physical_models/Reluctance.h"
#include "physical_models/InitialPermeability.h"
#include "physical_models/Impedance.h"
#include "physical_models/LeakageInductance.h"
#include "support/Settings.h"
#include <cmath>
#include <numbers>
#include <complex>
#include "support/Exceptions.h"


namespace OpenMagnetics {


Curve2D Sweeper::sweep_impedance_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, std::string mode, std::string title, bool fast) {
    std::vector<double> frequencies;
    if (mode == "linear") {
        frequencies = linear_spaced_array(start, stop, numberElements);
    }
    else if (mode == "log") {
        frequencies = logarithmic_spaced_array(start, stop, numberElements);
    }
    else {
        throw ModelNotAvailableException("Unknown spaced array mode");
    }

    // The terminal impedance is a series cascade of resonant tanks (a Foster
    // ladder): the magnetizing tank (first resonance) plus, for coupled magnetics,
    // one leakage tank per additional winding. The leakage inductance — the flux
    // that does not couple through the core — resonates with the inter-winding
    // capacitance to give the higher self-resonances (the common-mode-choke
    // "second resonance", generalized to any multi-winding magnetic). At low
    // frequency the leakage tanks reduce to ωL_leak ≪ magnetizing reactance, so a
    // single-winding inductor's single-resonance behaviour is preserved.
    //
    // The model holds the frequency-independent building blocks (reluctance, stray
    // capacitance, leakage, DC resistance) so they are computed ONCE here; only the
    // complex permeability, the (frequency-dependent) winding resistance and the
    // complex arithmetic run per point. With fast=true the leakage damping uses
    // DC + analytic skin effect; with fast=false it uses the full field-based
    // effective resistance (adds proximity), more accurate but slower per point.
    const double temperature = Defaults().ambientTemperature;
    double referenceFrequency = frequencies[frequencies.size() / 2];
    auto impedanceModel = OpenMagnetics::Impedance();
    auto model = impedanceModel.build_wideband_impedance_model(magnetic, referenceFrequency, temperature, fast);

    std::vector<double> impedances;
    impedances.reserve(frequencies.size());
    for (auto frequency : frequencies) {
        impedances.push_back(abs(impedanceModel.impedance_from_model(model, frequency)));
    }

    return Curve2D(frequencies, impedances, title);
}

Curve2D Sweeper::sweep_common_mode_impedance_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, std::string mode, std::string title) {
    std::vector<double> frequencies;
    if (mode == "linear") {
        frequencies = linear_spaced_array(start, stop, numberElements);
    }
    else if (mode == "log") {
        frequencies = logarithmic_spaced_array(start, stop, numberElements);
    }
    else {
        throw ModelNotAvailableException("Unknown spaced array mode");
    }

    // Common mode drives all windings in parallel, so the leakage between them is
    // not excited and the terminal impedance is the magnetizing tank alone: the
    // air-cored inductance times the core complex permeability µ(f) (µ'' is the
    // damping), shunted by the winding self-capacitance. The frequency-independent
    // building blocks (reluctance, stray capacitance) are computed ONCE here; only
    // µ(f) and the complex arithmetic run per point.
    auto impedanceModel = OpenMagnetics::Impedance();
    auto model = impedanceModel.build_common_mode_impedance_model(magnetic);

    std::vector<double> impedances;
    impedances.reserve(frequencies.size());
    for (auto frequency : frequencies) {
        impedances.push_back(abs(impedanceModel.impedance_from_model(model, frequency)));
    }

    return Curve2D(frequencies, impedances, title);
}

Curve2D Sweeper::sweep_differential_mode_impedance_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, std::string mode, std::string title) {
    std::vector<double> frequencies;
    if (mode == "linear") {
        frequencies = linear_spaced_array(start, stop, numberElements);
    }
    else if (mode == "log") {
        frequencies = logarithmic_spaced_array(start, stop, numberElements);
    }
    else {
        throw ModelNotAvailableException("Unknown spaced array mode");
    }

    // The leakage inductance, winding resistance and inter-winding capacitance are
    // frequency-independent, but the stray-capacitance model behind the last one is
    // expensive (O(N^2) over turns). Compute the parameters ONCE (leakage taken at a
    // mid-band reference) and evaluate the impedance cheaply at each point — otherwise
    // the whole capacitance model would run per frequency (~100x slower).
    double referenceFrequency = std::sqrt(frequencies.front() * frequencies.back());
    auto impedance = OpenMagnetics::Impedance();
    auto parameters = impedance.calculate_differential_mode_parameters(magnetic.get_core(), magnetic.get_coil(), referenceFrequency);

    std::vector<double> impedances;
    for (auto frequency : frequencies) {
        impedances.push_back(abs(impedance.differential_mode_impedance_from_parameters(parameters, frequency)));
    }

    return Curve2D(frequencies, impedances, title);
}

Curve2D Sweeper::sweep_q_factor_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, std::string mode, std::string title) {
    std::vector<double> frequencies;
    if (mode == "linear") {
        frequencies = linear_spaced_array(start, stop, numberElements);
    }
    else if (mode == "log") {
        frequencies = logarithmic_spaced_array(start, stop, numberElements);
    }
    else {
        throw ModelNotAvailableException("Unknown spaced array mode");
    }

    std::vector<double> qFactors;
    for (auto frequency : frequencies) {
        auto qFactor = abs(OpenMagnetics::Impedance().calculate_q_factor(magnetic, frequency));
        qFactors.push_back(qFactor);
    }

    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();
    return Curve2D(frequencies, qFactors, title);
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
        throw ModelNotAvailableException("Unknown spaced array mode");
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
        throw ModelNotAvailableException("Unknown spaced array mode");
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
        throw ModelNotAvailableException("Unknown spaced array mode");
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
        throw ModelNotAvailableException("Unknown spaced array mode");
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
        throw ModelNotAvailableException("Unknown spaced array mode");
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
        throw ModelNotAvailableException("Unknown spaced array mode");
    }
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();

    auto magnetizingInductanceModel = MagnetizingInductance();
    auto magnetizingInductance = resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(core, coil).get_magnetizing_inductance());

    std::vector<double> coreResistances;
    auto coreLossesModelSteinmetz = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Steinmetz"}}));
    auto coreLossesModelProprietary = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Proprietary"}}));
    auto coreLossesMethods = Core::get_available_core_losses_methods(core.resolve_material());

    if (std::find(coreLossesMethods.begin(), coreLossesMethods.end(), VolumetricCoreLossesMethodType::STEINMETZ) != coreLossesMethods.end()) {
        for (auto frequency : frequencies) {
            auto coreResistance =  coreLossesModelSteinmetz->get_core_losses_series_resistance(core, frequency, temperature, magnetizingInductance);
            coreResistances.push_back(coreResistance);
        }
    }
    else {
        for (auto frequency : frequencies) {
            auto coreResistance =  coreLossesModelProprietary->get_core_losses_series_resistance(core, frequency, temperature, magnetizingInductance);
            coreResistances.push_back(coreResistance);
        }
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
        throw ModelNotAvailableException("Unknown spaced array mode");
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
        auto voltageExcitation = Inputs::calculate_induced_voltage(excitation, magnetizingInductance);
        excitation.set_voltage(voltageExcitation);

        if (numberWindings == 1 && excitation.get_current()) {
            Inputs::set_current_as_magnetizing_current(&operatingPoint);
        }
        else if (Inputs::is_multiport_inductor(operatingPoint, coil.get_isolation_sides())) {
            auto magnetizingCurrent = Inputs::get_multiport_inductor_magnetizing_current(operatingPoint);
            excitation.set_magnetizing_current(magnetizingCurrent);
            operatingPoint.get_mutable_excitations_per_winding()[0] = excitation;
        }
        else if (Inputs::can_be_common_mode_choke(operatingPoint) && core.get_type() == CoreType::TOROIDAL) {
            // CMC: drive the core from the common-mode current alone, so the
            // DM (line) current which cancels in the toroid doesn't inflate B.
            // Same pattern MagnetizingInductance.cpp:149 uses.
            auto magnetizingCurrent = Inputs::get_common_mode_choke_magnetizing_current(operatingPoint);
            excitation.set_magnetizing_current(magnetizingCurrent);
            operatingPoint.get_mutable_excitations_per_winding()[0] = excitation;
        }
        else if (excitation.get_voltage()) {
            auto voltage = excitation.get_voltage().value();
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
        throw ModelNotAvailableException("Unknown spaced array mode");
    }

    auto magnetizingInductanceModel = MagnetizingInductance();
    auto turnsRatios = magnetic.get_mutable_coil().get_turns_ratios();
    auto magnetizingInductance = resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil()).get_magnetizing_inductance());

    std::vector<double> windingLossesPerFrequency;
    for (auto frequency : frequencies) {
        Inputs::scale_time_to_frequency(operatingPoint, frequency, true);
        operatingPoint = Inputs::process_operating_point(operatingPoint, magnetizingInductance);

        auto windingLosses =  WindingLosses().calculate_losses(magnetic, operatingPoint, temperature).get_winding_losses();

        windingLossesPerFrequency.push_back(windingLosses);
    }

    return Curve2D(frequencies, windingLossesPerFrequency, title);
}

std::vector<ScalarMatrixAtFrequency> Sweeper::sweep_resistance_matrix_over_frequency(
    Magnetic magnetic, 
    double start, 
    double stop, 
    size_t numberElements, 
    double temperature, 
    std::string mode) {
    
    std::vector<double> frequencies;
    if (mode == "linear") {
        frequencies = linear_spaced_array(start, stop, numberElements);
    }
    else if (mode == "log") {
        frequencies = logarithmic_spaced_array(start, stop, numberElements);
    }
    else {
        throw ModelNotAvailableException("Unknown spaced array mode");
    }

    std::vector<ScalarMatrixAtFrequency> resistanceMatrices;
    resistanceMatrices.reserve(numberElements);
    
    WindingLosses windingLossesModel;
    
    for (auto frequency : frequencies) {
        auto resistanceMatrix = windingLossesModel.calculate_resistance_matrix(magnetic, temperature, frequency);
        resistanceMatrices.push_back(resistanceMatrix);
    }

    return resistanceMatrices;
}

} // namespace OpenMagnetics
