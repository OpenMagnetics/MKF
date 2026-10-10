#include "processors/MagneticSimulator.h"
#include "physical_models/LeakageInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/Impedance.h"
#include "physical_models/ComplexPermeability.h"
#include "physical_models/Temperature.h"
#include "support/Settings.h"
#include <MAS.hpp>
#include "support/Exceptions.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <string>


namespace OpenMagnetics {

Mas MagneticSimulator::simulate(Mas mas, bool fastMode){
    return simulate(mas.get_mutable_inputs(), mas.get_mutable_magnetic(), fastMode);
}
Mas MagneticSimulator::simulate(const Inputs& inputs, const Magnetic& magnetic, bool fastMode){
    set_topology(inputs.get_design_requirements().get_topology());
    Mas mas;
    std::vector<Outputs> outputs;
    std::vector<OperatingPoint> simulatedOperatingPoints;
    mas.set_inputs(inputs);

    // Every winding must carry its own excitation (ABT #814). Deeper models
    // catch SOME of these — winding losses throw "Missing current excitation
    // for winding 2" once a later winding is unexcited — but a two-winding
    // magnetic given a single excitation used to run all the way through and
    // report complete, plausible losses for a component whose second winding
    // the caller never described. Silent under-specification is worse than a
    // refusal: the numbers look right and cannot be checked.
    {
        size_t numberWindings = magnetic.get_coil().get_functional_description().size();
        const auto& operatingPoints = mas.get_inputs().get_operating_points();
        for (size_t index = 0; index < operatingPoints.size(); ++index) {
            size_t numberExcitations = operatingPoints[index].get_excitations_per_winding().size();
            if (numberExcitations != numberWindings) {
                throw InvalidInputException(
                    ErrorCode::MISSING_DATA,
                    "Operating point " + std::to_string(index) + " has " +
                    std::to_string(numberExcitations) + " excitations for a magnetic with " +
                    std::to_string(numberWindings) + " windings; every winding needs exactly one excitation");
            }
        }
    }

    for (auto& operatingPoint : mas.get_mutable_inputs().get_mutable_operating_points()){
        Outputs output;
        InductanceOutput inductanceOutput;
        inductanceOutput.set_magnetizing_inductance(calculate_magnetizing_inductance(operatingPoint, magnetic));
        if (!fastMode) {
            if (magnetic.get_coil().get_functional_description().size() > 1) {
                inductanceOutput.set_leakage_inductance(calculate_leakage_inductance(operatingPoint, magnetic));
            }
        }
        output.set_inductance(inductanceOutput);
        // ABT #838: core and winding losses together, each at the temperature the COMBINED
        // dissipation produces (calculate_core_losses iterates both on one network). Fast mode
        // runs the same solve: the published core temperature is the same quantity in both.
        // The full network needs turn nodes; wind a local copy if the coil is bare (as
        // calculate_temperature does).
        Magnetic woundMagnetic = magnetic;
        if (!woundMagnetic.get_coil().get_turns_description()) {
            woundMagnetic.get_mutable_coil().wind();
        }
        std::optional<WindingLossesOutput> windingLossesOutput;
        std::optional<TemperatureOutput> temperatureOutput;
        output.set_core_losses(calculate_core_losses(operatingPoint, woundMagnetic, &windingLossesOutput, &temperatureOutput));
        output.set_winding_losses(windingLossesOutput.value());
        if (!fastMode) {
            // Full-network hot-spot from exactly the losses published above -- the same config
            // TemperatureConfig::fromSimulatedOutput builds, so the exported MAS carries the
            // temperature the UI's temperature map shows (ABT #906).
            output.set_temperature(temperatureOutput.value());
        }

        outputs.push_back(output);
        simulatedOperatingPoints.push_back(operatingPoint);
    }
    mas.get_mutable_inputs().set_operating_points(simulatedOperatingPoints);
    mas.set_magnetic(magnetic);
    mas.set_outputs(outputs);
    return mas;
}

MagnetizingInductanceOutput MagneticSimulator::calculate_magnetizing_inductance(OperatingPoint& operatingPoint, Magnetic magnetic){
    return _magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic, &operatingPoint);
}

LeakageInductanceOutput MagneticSimulator::calculate_leakage_inductance(OperatingPoint& operatingPoint, Magnetic magnetic){
    double frequency = operatingPoint.get_excitations_per_winding()[0].get_frequency();
    return calculate_leakage_inductance(magnetic, frequency);
}

LeakageInductanceOutput MagneticSimulator::calculate_leakage_inductance(Magnetic magnetic, double frequency){
    // Same shape as the public calculate_leakage_inductance WASM/Python API: one entry
    // per winding, index = winding index, 0 at the source (primary) slot, every value
    // referred to the primary. The previous hand-rolled loop emitted a secondaries-only
    // (N-1) array into the same MAS field, so consumers could not tell which shape they
    // were reading (web bug reports 125/139: the UI showed the primary's 0 as "Llk").
    return OpenMagnetics::LeakageInductance().calculate_leakage_inductance_all_windings(magnetic, frequency);
}

WindingLossesOutput MagneticSimulator::calculate_winding_losses(OperatingPoint& operatingPoint, Magnetic magnetic, std::optional<double> temperature){
    auto& settings = OpenMagnetics::Settings::GetInstance();
    // RAII: if calculate_losses throws, the manual restore below never ran and the
    // process kept mirroring=1 for every subsequent computation.
    SettingsGuard<int> mirroringGuard(settings, &Settings::get_magnetic_field_mirroring_dimension, &Settings::set_magnetic_field_mirroring_dimension, 1);

    double simulationTemperature;
    if (!temperature) {
        simulationTemperature = operatingPoint.get_conditions().get_ambient_temperature();
    }
    else {
        simulationTemperature = temperature.value();
    }
    WindingLosses windingLosses;
    // windingLosses.set_winding_losses_harmonic_amplitude_threshold(0.01);
    auto losses = windingLosses.calculate_losses(magnetic, operatingPoint, simulationTemperature);
    return losses;
}

TemperatureOutput MagneticSimulator::calculate_temperature(OperatingPoint& operatingPoint, Magnetic magnetic, const Outputs& output) {
    // The full thermal network needs turn nodes; wind a local copy if the coil is bare.
    if (!magnetic.get_coil().get_turns_description()) {
        magnetic.get_mutable_coil().wind();
    }
    auto config = TemperatureConfig::fromSimulatedOutput(operatingPoint, output);
    auto thermalResult = Temperature(magnetic, config).calculateTemperatures();

    TemperatureOutput temperatureOutput;
    temperatureOutput.set_initial_temperature(config.ambientTemperature);
    temperatureOutput.set_maximum_temperature(thermalResult.maximumTemperature);
    double totalLosses = config.coreLosses + config.windingLosses;
    if (totalLosses > 0) {
        temperatureOutput.set_bulk_thermal_resistance(
            (thermalResult.maximumTemperature - config.ambientTemperature) / totalLosses);
    }
    temperatureOutput.set_method_used("ThermalNetwork");
    temperatureOutput.set_origin(ResultOrigin::SIMULATION);
    return temperatureOutput;
}

CoreLossesOutput MagneticSimulator::calculate_core_losses(OperatingPoint& operatingPoint, Magnetic magnetic,
                                                          std::optional<WindingLossesOutput>* windingLossesOutput,
                                                          std::optional<TemperatureOutput>* temperatureOutput) {
    if (!operatingPoint.get_excitations_per_winding()[0].get_current()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Missing current in operating point");
    }
    // ABT #838: a wound coil heats the core, and the core's temperature sets the winding's: both
    // losses are evaluated at the temperatures their COMBINED dissipation produces, on the full
    // network (turn nodes carry the winding power). Without turns there are no winding losses to
    // compute and no turn nodes to put them in (the core advisers' sentinel-coil candidates): the
    // core is then solved on the core-only network with its own losses.
    const bool withWinding = magnetic.get_coil().get_turns_description().has_value();
    if (windingLossesOutput != nullptr && !withWinding) {
        throw CoilNotProcessedException("Winding losses asked of '" + magnetic.get_reference() +
                                        "', whose coil has no turns: wind it first");
    }

    const double ambientTemperature = operatingPoint.get_conditions().get_ambient_temperature();
    double coreTemperature = ambientTemperature;
    double windingTemperature = ambientTemperature;
    std::string history;
    double coreRelaxation = 1.0;
    double windingRelaxation = 1.0;
    std::optional<double> previousCoreResidual;
    std::optional<double> previousWindingResidual;
    const auto coreMaterial = magnetic.get_core().resolve_material();
    const auto curieTemperature = coreMaterial.get_curie_temperature();

    for (size_t iteration = 1; iteration <= kOperatingTemperatureMaximumIterations; ++iteration) {
        // Core losses at the core's temperature. The magnetizing model reads the core temperature
        // from the operating point's conditions (permeability vs T sets B for a current-driven
        // part), so it gets a COPY at T_core; the caller's operating point keeps its real ambient.
        OperatingPoint operatingPointAtCoreTemperature = operatingPoint;
        operatingPointAtCoreTemperature.get_mutable_conditions().set_ambient_temperature(coreTemperature);
        auto magneticFluxDensity = _magnetizingInductanceModel.calculate_inductance_and_magnetic_flux_density(
            magnetic, &operatingPointAtCoreTemperature).second;
        // Read the excitation after the inductance model: it stamps the magnetizing current it
        // derived (a CMC's common-mode current, a voltage's integral) into the operating point it
        // was given. A copy taken before it lacked that current whenever the operating point
        // arrived without one, and the loss-factor model then refused it ("Missing magnetizing current").
        OperatingPointExcitation excitation = operatingPointAtCoreTemperature.get_excitations_per_winding()[0];
        excitation.set_magnetic_flux_density(magneticFluxDensity);
        CoreLossesOutput coreLossesOutput = _coreLossesModel.calculate_core_losses(magnetic.get_core(), excitation, coreTemperature);
        coreLossesOutput.set_temperature(coreTemperature);

        // One network, every heat source that is modelled.
        std::optional<WindingLossesOutput> windingLosses;
        TemperatureConfig config;
        if (withWinding) {
            windingLosses = calculate_winding_losses(operatingPoint, magnetic, windingTemperature);
            Outputs losses;
            losses.set_core_losses(coreLossesOutput);
            losses.set_winding_losses(windingLosses.value());
            config = TemperatureConfig::fromSimulatedOutput(operatingPoint, losses);
        }
        else {
            config = TemperatureConfig::fromMasOperatingConditions(operatingPoint.get_conditions());
            config.plotSchematic = false;
            config.coreLosses = coreLossesOutput.get_core_losses();
            config.coreOnly = true;
        }
        auto thermalResult = Temperature(magnetic, config).calculateTemperatures();
        if (!thermalResult.maximumCoreTemperature) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT,
                "Thermal network of '" + magnetic.get_reference() + "' has no core node: no core temperature to "
                "evaluate the core losses at");
        }
        if (withWinding && !thermalResult.averageTurnTemperature) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT,
                "Thermal network of '" + magnetic.get_reference() + "' has no turn node: no winding temperature to "
                "evaluate the winding losses at");
        }
        const double nextCoreTemperature = thermalResult.maximumCoreTemperature.value();
        const double nextWindingTemperature = withWinding ? thermalResult.averageTurnTemperature.value() : ambientTemperature;
        history += " [" + std::to_string(iteration) + "] Pcore=" + std::to_string(coreLossesOutput.get_core_losses()) +
                   " W at " + std::to_string(coreTemperature) + " C";
        if (magneticFluxDensity.get_processed() && magneticFluxDensity.get_processed()->get_peak()) {
            history += " (Bpk " + std::to_string(magneticFluxDensity.get_processed()->get_peak().value()) + " T)";
        }
        if (withWinding) {
            history += ", Pcu=" + std::to_string(windingLosses->get_winding_losses()) + " W at " +
                       std::to_string(windingTemperature) + " C";
        }
        history += " -> core " + std::to_string(nextCoreTemperature) + " C";
        if (withWinding) {
            history += ", winding " + std::to_string(nextWindingTemperature) + " C";
        }
        history += ";";

        if (std::fabs(nextCoreTemperature - coreTemperature) < kOperatingTemperatureTolerance &&
            std::fabs(nextWindingTemperature - windingTemperature) < kOperatingTemperatureTolerance) {
            // Converged: the losses were evaluated at temperatures within the tolerance of the
            // ones these very losses produce. Publish exactly the set that was solved.
            operatingPoint.get_mutable_excitations_per_winding()[0] = excitation;
            if (windingLossesOutput != nullptr) {
                *windingLossesOutput = windingLosses;
            }
            if (temperatureOutput != nullptr) {
                TemperatureOutput temperature;
                temperature.set_initial_temperature(config.ambientTemperature);
                temperature.set_maximum_temperature(thermalResult.maximumTemperature);
                double totalLosses = config.coreLosses + config.windingLosses;
                if (totalLosses > 0) {
                    temperature.set_bulk_thermal_resistance(
                        (thermalResult.maximumTemperature - config.ambientTemperature) / totalLosses);
                }
                temperature.set_method_used("ThermalNetwork");
                temperature.set_origin(ResultOrigin::SIMULATION);
                *temperatureOutput = temperature;
            }
            return coreLossesOutput;
        }
        // Next iterate. Plain substitution (relaxation 1) is the fastest when the losses' feedback
        // on temperature is weak, but a ferrite whose losses FALL steeply with temperature makes
        // the map overshoot the fixed point and swing around it forever (seen: 3.9 W at 106 C,
        // 2.6 W at 134 C, 0.01 W at 150 C). Each time a temperature's residual changes sign, that
        // temperature's step is halved, which turns the swing into convergence. Core and winding
        // keep their own step: damping the winding because the CORE swung would leave it crawling
        // toward its own answer. The convergence test above is on the residuals themselves, so
        // damping never makes an unsettled state look settled.
        const double coreResidual = nextCoreTemperature - coreTemperature;
        const double windingResidual = nextWindingTemperature - windingTemperature;
        if (previousCoreResidual && previousCoreResidual.value() * coreResidual < 0) {
            coreRelaxation *= 0.5;
        }
        if (previousWindingResidual && previousWindingResidual.value() * windingResidual < 0) {
            windingRelaxation *= 0.5;
        }
        previousCoreResidual = coreResidual;
        previousWindingResidual = windingResidual;
        double proposedCoreTemperature = coreTemperature + coreRelaxation * coreResidual;
        double proposedWindingTemperature = windingTemperature + windingRelaxation * windingResidual;

        // The core-loss and permeability models describe the FERROMAGNETIC material and refuse a
        // temperature at or above its Curie point. The first steps start from ambient-temperature
        // losses and can overshoot far past where the answer lies, so the iterate is kept inside
        // the range where the losses exist: at most halfway from the current core temperature to
        // the Curie point. If the iteration has been driven to within the tolerance of the Curie
        // point and the losses there still push the core hotter, no operating point was found
        // below it: refuse, saying so.
        if (curieTemperature && proposedCoreTemperature >= curieTemperature.value()) {
            const double distanceToCurie = curieTemperature.value() - coreTemperature;
            if (distanceToCurie < kOperatingTemperatureTolerance) {
                throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT,
                    "Losses of '" + magnetic.get_reference() + "': with the core at " + std::to_string(coreTemperature) +
                    " C, just below the Curie temperature of " + coreMaterial.get_name() + " (" +
                    std::to_string(curieTemperature.value()) + " C), its core and winding losses still drive it to " +
                    std::to_string(nextCoreTemperature) + " C. No operating temperature was found where the core is "
                    "ferromagnetic. History:" + history);
            }
            proposedCoreTemperature = coreTemperature + 0.5 * distanceToCurie;
        }
        coreTemperature = proposedCoreTemperature;
        windingTemperature = proposedWindingTemperature;
    }
    throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT,
        "Losses of '" + magnetic.get_reference() + "' and the temperatures they produce did not settle within " +
        std::to_string(kOperatingTemperatureMaximumIterations) + " iterations (tolerance " +
        std::to_string(kOperatingTemperatureTolerance) + " K): no self-consistent operating temperature exists "
        "(thermal runaway?) or the iteration oscillates. History:" + history);
}

namespace {
    // Builds a single-valued DimensionWithTolerance (nominal only) tagged with a unit.
    DimensionWithTolerance nominal_dimension(double value, const std::string& unit) {
        DimensionWithTolerance dimension;
        dimension.set_nominal(value);
        dimension.set_unit(unit);
        return dimension;
    }

    // RMS of a winding's current/voltage signal in an operating point, if present.
    std::optional<double> excitation_rms(const SignalDescriptor& signal) {
        if (!signal.get_processed() || !signal.get_processed()->get_rms()) {
            return std::nullopt;
        }
        return signal.get_processed()->get_rms();
    }
}

MagneticManufacturerInfo MagneticSimulator::build_datasheet(Mas& mas) {
    auto& magnetic = mas.get_mutable_magnetic();
    const auto& outputs = mas.get_outputs();
    const auto& operatingPoints = mas.get_inputs().get_operating_points();

    if (outputs.empty()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA,
            "Cannot build datasheet: MAS carries no simulation outputs; run MagneticSimulator::simulate first");
    }
    if (operatingPoints.size() != outputs.size()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA,
            "Cannot build datasheet: number of operating points (" + std::to_string(operatingPoints.size()) +
            ") does not match number of outputs (" + std::to_string(outputs.size()) + ")");
    }

    const auto& windings = magnetic.get_coil().get_functional_description();
    size_t numberWindings = windings.size();

    // ----------------------------------------------------------------- Part
    Part part;
    part.set_number_of_windings(static_cast<int64_t>(numberWindings));
    auto materialName = magnetic.get_mutable_core().get_material_name();
    if (!materialName.empty()) {
        part.set_material(materialName);
    }
    if (magnetic.get_manufacturer_info()) {
        auto existing = magnetic.get_manufacturer_info().value();
        if (existing.get_reference())   part.set_part_number(existing.get_reference());
        if (existing.get_description()) part.set_description(existing.get_description());
    }

    // ----------------------------------------------------------- Electrical
    MagneticDatasheetElectrical electrical;
    electrical.set_subtype(numberWindings > 1 ? ElectricalSubtype::COUPLED_INDUCTOR : ElectricalSubtype::INDUCTOR);

    // Magnetizing inductance and leakage inductance are taken from the nominal (first) operating point.
    if (outputs[0].get_inductance()) {
        auto inductanceOutput = outputs[0].get_inductance().value();
        electrical.set_inductance(inductanceOutput.get_magnetizing_inductance().get_magnetizing_inductance());

        // Per-winding leakage array is winding-indexed with 0 at the primary slot, so the
        // first secondary lives at index 1.
        if (inductanceOutput.get_leakage_inductance() &&
            inductanceOutput.get_leakage_inductance()->get_leakage_inductance_per_winding().size() > 1) {
            electrical.set_leakage_inductance(
                inductanceOutput.get_leakage_inductance()->get_leakage_inductance_per_winding()[1]);
        }
    }

    // DC resistance per winding, evaluated at the nominal operating point's ambient temperature, via MKF.
    double dcResistanceTemperature = operatingPoints[0].get_conditions().get_ambient_temperature();
    auto dcResistancePerWinding = WindingOhmicLosses::calculate_dc_resistance_per_winding(
        magnetic.get_coil(), dcResistanceTemperature);
    if (numberWindings == 1 && dcResistancePerWinding.size() == 1) {
        electrical.set_dc_resistance(nominal_dimension(dcResistancePerWinding[0], "Ohm"));
    }
    else if (dcResistancePerWinding.size() == numberWindings) {
        std::vector<DimensionWithTolerance> dcResistances;
        for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
            dcResistances.push_back(nominal_dimension(dcResistancePerWinding[windingIndex], "Ohm"));
        }
        electrical.set_dc_resistances(dcResistances);
    }

    // Rated current: the heat-limited DC current for the standard 40 K rise (datasheet Irms),
    // computed from the winding ohmic loss and the core temperature model.
    electrical.set_rated_currents(std::vector<double>{magnetic.calculate_rated_current()});

    // Self-resonant frequency and an impedance-vs-frequency sweep. The sweep uses the wideband
    // Foster-ladder model (the magnetizing resonance plus one leakage resonance per extra
    // winding), built once and evaluated per frequency; stray capacitance and winding
    // resistance are folded in by the model. Anchored on the SRF and spanning several decades
    // above it so the published curve captures the higher (leakage) resonances too.
    {
        Impedance impedanceModel;
        double selfResonantFrequency = impedanceModel.calculate_self_resonant_frequency(magnetic, dcResistanceTemperature);
        electrical.set_self_resonant_frequency(selfResonantFrequency);

        auto widebandModel = impedanceModel.build_wideband_impedance_model(magnetic, selfResonantFrequency, dcResistanceTemperature);

        const size_t numberPoints = 51;
        double startFrequency = selfResonantFrequency / 1000;
        double stopFrequency = selfResonantFrequency * 100;
        // The published curve covers only the frequencies where the core material has complex
        // permeability data: outside its table there is nothing to evaluate (get_complex_permeability
        // throws there). This sweep is the datasheet's own choice of band, so it asks for the
        // material's range explicitly and keeps the part of the band that lies inside it. If none
        // of the band does, the curve cannot be computed and that is reported, not filled in.
        {
            auto coreMaterial = magnetic.get_core().resolve_material();
            auto [minimumMaterialFrequency, maximumMaterialFrequency] = ComplexPermeability().get_frequency_range(coreMaterial);
            double requestedStartFrequency = startFrequency;
            double requestedStopFrequency = stopFrequency;
            startFrequency = std::max(startFrequency, minimumMaterialFrequency);
            stopFrequency = std::min(stopFrequency, maximumMaterialFrequency);
            if (!(startFrequency < stopFrequency)) {
                throw ComplexPermeabilityFrequencyOutOfRangeException(coreMaterial.get_name(),
                    requestedStartFrequency < minimumMaterialFrequency ? requestedStopFrequency : requestedStartFrequency,
                    minimumMaterialFrequency, maximumMaterialFrequency);
            }
        }
        // Exact ends (logarithmic_spaced_array): start * 10^(n * step) can land one ulp past the
        // material's last tabulated frequency, where complex permeability throws.
        auto frequencies = logarithmic_spaced_array(startFrequency, stopFrequency, numberPoints);

        std::vector<DatasheetImpedancePoint> impedancePoints;
        double maximumImpedanceMagnitude = std::numeric_limits<double>::lowest();
        for (size_t pointIndex = 0; pointIndex < numberPoints; ++pointIndex) {
            double frequency = frequencies[pointIndex];
            auto impedance = impedanceModel.impedance_from_model(widebandModel, frequency);

            ImpedancePoint impedancePoint;
            impedancePoint.set_magnitude(std::abs(impedance));
            impedancePoint.set_real_part(impedance.real());
            impedancePoint.set_imaginary_part(impedance.imag());
            impedancePoint.set_phase(std::arg(impedance));

            DatasheetImpedancePoint datasheetImpedancePoint;
            datasheetImpedancePoint.set_frequency(frequency);
            datasheetImpedancePoint.set_impedance(impedancePoint);
            impedancePoints.push_back(datasheetImpedancePoint);

            maximumImpedanceMagnitude = std::max(maximumImpedanceMagnitude, std::abs(impedance));
        }
        electrical.set_impedance_points(impedancePoints);
        electrical.set_maximum_impedance(maximumImpedanceMagnitude);
    }

    // Turns ratio (primary : first secondary) for coupled parts.
    if (numberWindings > 1) {
        auto turnsRatios = magnetic.get_turns_ratios();
        if (!turnsRatios.empty()) {
            electrical.set_turns_ratios(std::vector<DimensionWithTolerance>{nominal_dimension(turnsRatios[0], "")});
        }
    }

    // Peak saturation current, evaluated at the hottest operating point (lowest, conservative I_sat),
    // delegating to the MKF saturation-current model.
    size_t hottestOperatingPointIndex = 0;
    double hottestAmbient = std::numeric_limits<double>::lowest();
    for (size_t operatingPointIndex = 0; operatingPointIndex < operatingPoints.size(); ++operatingPointIndex) {
        double ambient = operatingPoints[operatingPointIndex].get_conditions().get_ambient_temperature();
        if (ambient > hottestAmbient) {
            hottestAmbient = ambient;
            hottestOperatingPointIndex = operatingPointIndex;
        }
    }
    {
        OperatingPoint saturationOperatingPoint = operatingPoints[hottestOperatingPointIndex];
        electrical.set_saturation_current_peak(
            magnetic.calculate_saturation_current(saturationOperatingPoint, hottestAmbient));
    }

    // -------------------------------------------------------------- Thermal
    Thermal thermal;

    // Operating temperature range = span of ambient conditions across operating points.
    double minimumAmbient = std::numeric_limits<double>::max();
    double maximumAmbient = std::numeric_limits<double>::lowest();
    for (const auto& operatingPoint : operatingPoints) {
        double ambient = operatingPoint.get_conditions().get_ambient_temperature();
        minimumAmbient = std::min(minimumAmbient, ambient);
        maximumAmbient = std::max(maximumAmbient, ambient);
    }
    {
        DimensionWithTolerance operatingTemperature;
        operatingTemperature.set_minimum(minimumAmbient);
        operatingTemperature.set_maximum(maximumAmbient);
        operatingTemperature.set_nominal(operatingPoints[0].get_conditions().get_ambient_temperature());
        operatingTemperature.set_unit("Celsius");
        thermal.set_operating_temperature(operatingTemperature);
    }

    // Worst-case hot-spot temperature rise and thermal resistance, computed with the full
    // Temperature thermal-network model from the simulated core + winding losses at each
    // operating point (the per-turn winding-loss distribution is passed through for fidelity).
    Magnetic thermalMagnetic = magnetic;
    if (!thermalMagnetic.get_coil().get_turns_description()) {
        thermalMagnetic.get_mutable_coil().wind();
    }
    double worstTemperatureRise = std::numeric_limits<double>::lowest();
    std::optional<double> worstThermalResistance;
    for (size_t operatingPointIndex = 0; operatingPointIndex < outputs.size(); ++operatingPointIndex) {
        if (!outputs[operatingPointIndex].get_core_losses()) {
            continue;
        }
        double ambient = operatingPoints[operatingPointIndex].get_conditions().get_ambient_temperature();
        // Shared config builder (ABT #906) — same as outputs[].temperature and the plot wrappers.
        auto config = TemperatureConfig::fromSimulatedOutput(operatingPoints[operatingPointIndex],
                                                             outputs[operatingPointIndex]);

        double hotspot = Temperature(thermalMagnetic, config).calculateTemperatures().maximumTemperature;
        double temperatureRise = hotspot - ambient;
        if (temperatureRise > worstTemperatureRise) {
            worstTemperatureRise = temperatureRise;
            double totalLosses = config.coreLosses + config.windingLosses;
            worstThermalResistance = totalLosses > 0 ? std::optional<double>(temperatureRise / totalLosses)
                                                     : std::nullopt;
        }
    }
    if (worstTemperatureRise > std::numeric_limits<double>::lowest()) {
        thermal.set_temperature_rise(worstTemperatureRise);
        if (worstThermalResistance) {
            thermal.set_thermal_resistance(worstThermalResistance.value());
        }
    }

    // ----------------------------------------------------------- Mechanical
    Mechanical mechanical;
    auto dimensions = magnetic.get_maximum_dimensions();  // [width, height, depth] in metres
    mechanical.set_width(nominal_dimension(dimensions[0], "m"));
    mechanical.set_height(nominal_dimension(dimensions[1], "m"));
    mechanical.set_length(nominal_dimension(dimensions[2], "m"));

    // ---------------------------------------------------------- Application
    // The application block states the circuit conditions this design was simulated for.
    MagneticDatasheetApplication application;
    application.set_switching_frequency(operatingPoints[0].get_excitations_per_winding()[0].get_frequency());

    // Input voltage range (primary winding) across operating points.
    double minimumInputVoltage = std::numeric_limits<double>::max();
    double maximumInputVoltage = std::numeric_limits<double>::lowest();
    for (const auto& operatingPoint : operatingPoints) {
        const auto& excitation = operatingPoint.get_excitations_per_winding()[0];
        if (excitation.get_voltage()) {
            auto rms = excitation_rms(excitation.get_voltage().value());
            if (rms) {
                minimumInputVoltage = std::min(minimumInputVoltage, rms.value());
                maximumInputVoltage = std::max(maximumInputVoltage, rms.value());
            }
        }
    }
    if (maximumInputVoltage > std::numeric_limits<double>::lowest()) {
        DimensionWithTolerance inputVoltage;
        inputVoltage.set_minimum(minimumInputVoltage);
        inputVoltage.set_maximum(maximumInputVoltage);
        inputVoltage.set_unit("V");
        application.set_input_voltage(Dimension(inputVoltage));
    }

    // Output voltages (secondary windings) at the nominal operating point.
    if (numberWindings > 1 &&
        operatingPoints[0].get_excitations_per_winding().size() == numberWindings) {
        std::vector<double> outputVoltages;
        for (size_t windingIndex = 1; windingIndex < numberWindings; ++windingIndex) {
            const auto& excitation = operatingPoints[0].get_excitations_per_winding()[windingIndex];
            if (excitation.get_voltage()) {
                auto rms = excitation_rms(excitation.get_voltage().value());
                if (rms) {
                    outputVoltages.push_back(rms.value());
                }
            }
        }
        if (!outputVoltages.empty()) {
            application.set_output_voltages(outputVoltages);
        }
    }

    // ------------------------------------------------- Assemble datasheetInfo
    DatasheetInfo datasheetInfo;
    datasheetInfo.set_part(part);
    datasheetInfo.set_electrical(std::vector<MagneticDatasheetElectrical>{electrical});
    datasheetInfo.set_thermal(thermal);
    datasheetInfo.set_mechanical(mechanical);
    datasheetInfo.set_application(application);

    // Preserve any existing manufacturer metadata; only (re)populate the datasheetInfo block.
    MagneticManufacturerInfo manufacturerInfo;
    if (magnetic.get_manufacturer_info()) {
        manufacturerInfo = magnetic.get_manufacturer_info().value();
    }
    else {
        manufacturerInfo.set_name("OpenMagnetics");
    }
    manufacturerInfo.set_datasheet_info(datasheetInfo);

    magnetic.set_manufacturer_info(manufacturerInfo);
    return manufacturerInfo;
}

} // namespace OpenMagnetics
