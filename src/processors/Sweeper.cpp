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
#include <algorithm>
#include <limits>
#include <optional>
#include <sstream>
#include <numbers>
#include <complex>
#include "support/Exceptions.h"


namespace OpenMagnetics {


Curve2D Sweeper::sweep_impedance_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, std::string mode, std::string title, bool fast, bool fastCapacitance) {
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
    auto impedanceModel = OpenMagnetics::Impedance(fastCapacitance);
    auto model = impedanceModel.build_wideband_impedance_model(magnetic, referenceFrequency, temperature, fast);

    std::vector<double> impedances;
    impedances.reserve(frequencies.size());
    for (auto frequency : frequencies) {
        impedances.push_back(abs(impedanceModel.impedance_from_model(model, frequency)));
    }

    return Curve2D(frequencies, impedances, title);
}

Curve2D Sweeper::sweep_common_mode_impedance_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, std::string mode, std::string title, bool fastCapacitance) {
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
    auto impedanceModel = OpenMagnetics::Impedance(fastCapacitance);
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
    auto parameters = impedance.calculate_differential_mode_parameters(magnetic, referenceFrequency);

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

    // Saturated-state model of this sweep (ABT #1790, Alf 2026-10-12). MagnetizingInductance refuses a
    // bias that saturates the core (CALCULATION_CORE_SATURATES, ABT #1740), and everywhere else that
    // refusal stands. A DC-bias sweep is drawn to see where the inductance collapses, so here the curve
    // continues past that point with the core treated as air: mu_r = 1 in the same reluctance model
    // (ungapped core path at mu_r = 1 in series with the gaps), L = N^2 / R. Those points are flagged in
    // the curve and named in its title. Only that error code is taken; anything else propagates.
    std::optional<double> airCoreMagnetizingInductance;
    auto get_air_core_magnetizing_inductance = [&]() {
        if (!airCoreMagnetizingInductance) {
            double numberTurns = static_cast<double>(inductorCoil.get_functional_description()[0].get_number_turns());
            double airCoreReluctance = ReluctanceModel::factory()->get_core_reluctance(magnetic.get_core(), 1.0).get_core_reluctance();
            airCoreMagnetizingInductance = numberTurns * numberTurns / airCoreReluctance;
        }
        return airCoreMagnetizingInductance.value();
    };

    std::vector<double> magnetizingInductances;
    std::vector<bool> saturatedPoints;
    std::optional<double> firstSaturatedCurrentOffset;
    for (auto currentOffset : currentOffsets) {
        auto operatingPoint = Inputs::create_operating_point_with_sinusoidal_current_mask(defaults.measurementFrequency, staticMagnetizingInductance, temperature, turnsRatios, currentMask, currentOffset);
        double magnetizingInductance;
        bool saturated = false;
        try {
            magnetizingInductance = resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), inductorCoil, &operatingPoint).get_magnetizing_inductance());
        }
        catch (const OpenMagneticsException& exception) {
            if (exception.code() != ErrorCode::CALCULATION_CORE_SATURATES) {
                throw;
            }
            saturated = true;
            magnetizingInductance = get_air_core_magnetizing_inductance();
            if (!firstSaturatedCurrentOffset) {
                firstSaturatedCurrentOffset = currentOffset;
            }
        }
        magnetizingInductances.push_back(magnetizingInductance);
        saturatedPoints.push_back(saturated);
    }

    if (firstSaturatedCurrentOffset) {
        std::ostringstream note;
        note << " (core saturated from " << firstSaturatedCurrentOffset.value()
             << " A: those points are the air-core inductance, mu_r = 1, "
             << get_air_core_magnetizing_inductance() << " H)";
        title += note.str();
    }

    Curve2D curve(currentOffsets, magnetizingInductances, title);
    curve.set_saturated_points(saturatedPoints);
    return curve;
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

// The loss model sweep_core_resistance_over_frequency evaluates for a material.
static std::string core_resistance_losses_model_name(const CoreMaterial& material) {
    // ABT #388: the model was chosen as "Steinmetz if the material has Steinmetz, else
    // Proprietary" — a two-way choice that silently assumed every non-Steinmetz material is a
    // vendor-proprietary one. A LOSS-FACTOR material (tan_delta/mu_i vs frequency, which is how
    // NiZn makers publish loss and often the ONLY thing they publish) therefore landed in the
    // proprietary branch and was rejected with "No proprietary volumetric losses method",
    // blocking the whole subcircuit export for it — 101 of 118 shielded-drum models on just two
    // datasheet-sourced grades. MKF already has a LossFactor model; the selection simply never
    // offered it. Ask for the model that fits the material instead of guessing between two.
    auto coreLossesMethods = Core::get_available_core_losses_methods(material);
    auto hasMethod = [&](VolumetricCoreLossesMethodType method) {
        return std::find(coreLossesMethods.begin(), coreLossesMethods.end(), method) != coreLossesMethods.end();
    };
    if (hasMethod(VolumetricCoreLossesMethodType::STEINMETZ)) {
        return "Steinmetz";
    }
    if (hasMethod(VolumetricCoreLossesMethodType::LOSS_FACTOR)) {
        return "LossFactor";
    }
    return "Proprietary";
}

std::pair<double, double> Sweeper::core_resistance_frequency_window(Magnetic magnetic, double start, double stop) {
    auto material = magnetic.get_core().resolve_material();
    if (core_resistance_losses_model_name(material) != "Steinmetz") {
        return {start, stop};
    }
    auto [spanMinimum, spanMaximum] = CoreLossesModel::get_steinmetz_fitted_span(material);
    // Settings::allowMaterialDataExtrapolation (explicit opt-in): sweep the whole requested window; the
    // points outside the span are extrapolated and warned about by get_steinmetz_coefficients.
    if (Settings::GetInstance().get_allow_material_data_extrapolation()) {
        return {start, stop};
    }
    double clippedStart = std::max(start, spanMinimum);
    double clippedStop = std::min(stop, spanMaximum);
    if (!(clippedStop > clippedStart)) {
        throw MaterialFrequencyOutOfSpanException(material.get_name(), "Steinmetz",
                                                  stop < spanMinimum ? stop : start, spanMinimum, spanMaximum);
    }
    return {clippedStart, clippedStop};
}

Curve2D Sweeper::sweep_core_resistance_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, double temperature, std::string mode, std::string title) {
    // ABT #1456: a Steinmetz material is only evaluated inside its fitted span; outside it the coefficients
    // throw. Sweep the part of [start, stop] the fit covers and say so in the title.
    auto [windowStart, windowStop] = core_resistance_frequency_window(magnetic, start, stop);
    if (windowStart != start || windowStop != stop) {
        title += " (restricted to the material's fitted Steinmetz span: " + std::to_string(windowStart) + " Hz to " +
                 std::to_string(windowStop) + " Hz)";
        start = windowStart;
        stop = windowStop;
    }
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
    // The window ends are the fitted span's ends: pin them exactly, since exp(log(f)) can land an ulp outside
    // the span and throw (ABT #1456).
    if (frequencies.empty()) {
        throw InvalidInputException("A frequency sweep needs at least 1 point");
    }
    frequencies.front() = start;
    if (frequencies.size() > 1) {
        frequencies.back() = stop;
    }
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();

    auto magnetizingInductanceModel = MagnetizingInductance();
    auto magnetizingInductance = resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(core, coil).get_magnetizing_inductance());

    std::vector<double> coreResistances;
    std::string coreLossesModelName = core_resistance_losses_model_name(core.resolve_material());
    auto coreLossesModel = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", coreLossesModelName}}));

    for (auto frequency : frequencies) {
        auto coreResistance = coreLossesModel->get_core_losses_series_resistance(core, frequency, temperature, magnetizingInductance);
        coreResistances.push_back(coreResistance);
    }


    return Curve2D(frequencies, coreResistances, title);
}

Curve2D Sweeper::sweep_core_losses_over_frequency(Magnetic magnetic, OperatingPoint operatingPoint, std::optional<Topology> topology, double start, double stop, size_t numberElements, double temperature, std::string mode, std::string title) {
    // ABT #1456: this sweep evaluates the Steinmetz model, whose coefficients exist only inside the material's
    // fitted span. Sweep the part of [start, stop] the fit covers and say so in the title, as the core
    // resistance sweep does; a request entirely outside the span throws.
    auto [windowStart, windowStop] = core_resistance_frequency_window(magnetic, start, stop);
    if (windowStart != start || windowStop != stop) {
        title += " (restricted to the material's fitted Steinmetz span: " + std::to_string(windowStart) + " Hz to " +
                 std::to_string(windowStop) + " Hz)";
        start = windowStart;
        stop = windowStop;
    }
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
    // The window ends are the fitted span's ends: pin them exactly, since exp(log(f)) can land an ulp outside
    // the span and throw (ABT #1456).
    if (frequencies.empty()) {
        throw InvalidInputException("A frequency sweep needs at least 1 point");
    }
    frequencies.front() = start;
    if (frequencies.size() > 1) {
        frequencies.back() = stop;
    }
    auto core = magnetic.get_core();

    // Each point's magnetizing current and flux density come from MagnetizingInductance, the model
    // the simulator and the advisers use, given the design's topology. This sweep kept its own copy
    // of that model's branches, which drifted from it: once the model told a multi-winding DMC from
    // a CMC, the copy still drove a DMC with the CMC triangle of its line current's whole swing at
    // the ripple frequency (ABT #1754).
    MagnetizingInductance magnetizingInductanceModel;
    magnetizingInductanceModel.set_topology(topology);
    operatingPoint.get_mutable_conditions().set_ambient_temperature(temperature);

    std::vector<double> coreLossesPerFrequency;
    // Saturated-state model of this sweep (ABT #1790, Alf 2026-10-12). The sweep holds the applied
    // voltage, so the flux swing grows as the frequency falls and the low end of a sweep can ask for
    // volt-seconds the core cannot take; MagnetizingInductance refuses those (CALCULATION_CORE_SATURATES,
    // ABT #1740). A saturated core is treated as air, and no core loss is computed for air or read from a
    // loss model past saturation: such a point stays in the curve at its frequency, flagged, with no
    // value (NaN), and the title names the saturated points. Only that error code is taken.
    std::vector<bool> saturatedPoints;
    std::vector<double> saturatedFrequencies;
    CoreLosses coreLossesModel;
    coreLossesModel.set_core_losses_model_name(CoreLossesModels::STEINMETZ);

    for (auto frequency : frequencies) {
        // Keep BOTH signals through the rescale. The third argument used to be `true`
        // (cleanFrequencyDependentFields), which with the default useCurrentAsBase=true
        // DELETES the voltage and keeps the current — so every point on the curve was
        // evaluated at the same magnetizing current, hence the same flux density, and
        // the losses climbed as f^alpha. A real design is voltage-driven: raise the
        // switching frequency at a fixed applied voltage and B falls as 1/f, so its
        // core losses go DOWN, not up. The two only ever agreed at the operating
        // point's own frequency, which is why the graph read tens of watts where the
        // Core Info panel and the datasheet read a fraction of one (ABT: user report).
        auto sweptOperatingPoint = operatingPoint;
        // Recompute the signals' harmonics and processed data on the rescaled time axis: the stored
        // ones describe the operating point's own frequency (ABT #1764).
        Inputs::scale_time_to_frequency(sweptOperatingPoint, frequency, false, true);
        // The model derives the magnetizing current and flux density at each frequency, as this
        // sweep always did: one the operating point carries may come from a model that did not
        // know the topology (the CMC triangle a DMC used to get).
        for (auto& excitation : sweptOperatingPoint.get_mutable_excitations_per_winding()) {
            excitation.set_magnetizing_current(std::nullopt);
            excitation.set_magnetic_flux_density(std::nullopt);
        }

        std::optional<SignalDescriptor> magneticFluxDensity;
        try {
            magneticFluxDensity = magnetizingInductanceModel.calculate_inductance_and_magnetic_flux_density(magnetic, &sweptOperatingPoint).second;
        }
        catch (const OpenMagneticsException& exception) {
            if (exception.code() != ErrorCode::CALCULATION_CORE_SATURATES) {
                throw;
            }
        }
        if (!magneticFluxDensity) {
            coreLossesPerFrequency.push_back(std::numeric_limits<double>::quiet_NaN());
            saturatedPoints.push_back(true);
            saturatedFrequencies.push_back(frequency);
            continue;
        }
        auto excitation = sweptOperatingPoint.get_excitations_per_winding()[0];
        excitation.set_magnetic_flux_density(magneticFluxDensity.value());

        auto coreLosses = coreLossesModel.calculate_core_losses(core, excitation, temperature).get_core_losses();
        coreLossesPerFrequency.push_back(coreLosses);
        saturatedPoints.push_back(false);
    }

    if (!saturatedFrequencies.empty()) {
        std::ostringstream note;
        note << " (core saturated at " << saturatedFrequencies.size() << " of " << frequencies.size() << " points, "
             << *std::min_element(saturatedFrequencies.begin(), saturatedFrequencies.end()) << " Hz to "
             << *std::max_element(saturatedFrequencies.begin(), saturatedFrequencies.end())
             << " Hz: the core is air there (mu_r = 1) and no core loss is computed)";
        title += note.str();
    }

    Curve2D curve(frequencies, coreLossesPerFrequency, title);
    curve.set_saturated_points(saturatedPoints);
    return curve;
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
        // Rescale the time axis and recompute the harmonics the loss model needs, but
        // do NOT reshape the excitation. process_operating_point() reflects waveforms,
        // re-derives the magnetizing current and synthesises a voltage; fed an operating
        // point whose voltage had just been deleted by the rescale, it handed the loss
        // model different currents from the ones the Winding Losses panel uses, and the
        // curve sat ~44% above the panel even at the operating point's own frequency.
        Inputs::scale_time_to_frequency(operatingPoint, frequency, false, true);

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
