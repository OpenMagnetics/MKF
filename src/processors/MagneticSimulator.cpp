#include "processors/MagneticSimulator.h"
#include "physical_models/LeakageInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/Impedance.h"
#include "physical_models/Temperature.h"
#include "support/Settings.h"
#include <MAS.hpp>
#include "support/Exceptions.h"

#include <cmath>
#include <complex>
#include <limits>
#include <string>


namespace OpenMagnetics {

Mas MagneticSimulator::simulate(Mas mas, bool fastMode){
    return simulate(mas.get_mutable_inputs(), mas.get_mutable_magnetic(), fastMode);
}
Mas MagneticSimulator::simulate(const Inputs& inputs, const Magnetic& magnetic, bool fastMode){
    Mas mas;
    std::vector<Outputs> outputs;
    std::vector<OperatingPoint> simulatedOperatingPoints;
    mas.set_inputs(inputs);

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
        output.set_core_losses(calculate_core_losses(operatingPoint, magnetic));
        output.set_winding_losses(calculate_winding_losses(operatingPoint, magnetic, operatingPoint.get_conditions().get_ambient_temperature()));

        outputs.push_back(output);
        simulatedOperatingPoints.push_back(operatingPoint);
    }
    mas.get_mutable_inputs().set_operating_points(simulatedOperatingPoints);
    mas.set_magnetic(magnetic);
    mas.set_outputs(outputs);
    return mas;
}

MagnetizingInductanceOutput MagneticSimulator::calculate_magnetizing_inductance(OperatingPoint& operatingPoint, Magnetic magnetic){
    return _magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil(), &operatingPoint);
}

LeakageInductanceOutput MagneticSimulator::calculate_leakage_inductance(OperatingPoint& operatingPoint, Magnetic magnetic){
    double frequency = operatingPoint.get_excitations_per_winding()[0].get_frequency();
    return calculate_leakage_inductance(magnetic, frequency);
}

LeakageInductanceOutput MagneticSimulator::calculate_leakage_inductance(Magnetic magnetic, double frequency){
    LeakageInductanceOutput leakageInductanceOutput; 
    for (size_t windingIndex = 1; windingIndex < magnetic.get_coil().get_functional_description().size(); ++windingIndex) {
        auto aux = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 0, windingIndex);
        if (windingIndex == 1) {
            leakageInductanceOutput = aux;
            leakageInductanceOutput.set_leakage_inductance_per_winding(std::vector<DimensionWithTolerance>());
        }
        auto currentLeakageInductancePerWinding = leakageInductanceOutput.get_leakage_inductance_per_winding();
        currentLeakageInductancePerWinding.push_back(aux.get_leakage_inductance_per_winding()[0]);
        leakageInductanceOutput.set_leakage_inductance_per_winding(currentLeakageInductancePerWinding);
    }
    return leakageInductanceOutput;
}

WindingLossesOutput MagneticSimulator::calculate_winding_losses(OperatingPoint& operatingPoint, Magnetic magnetic, std::optional<double> temperature){
    auto& settings = OpenMagnetics::Settings::GetInstance();
    auto oldMagneticFieldMirroringDimension = settings.get_magnetic_field_mirroring_dimension();
    settings.set_magnetic_field_mirroring_dimension(1);

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
    settings.set_magnetic_field_mirroring_dimension(oldMagneticFieldMirroringDimension);
    return losses;
}

CoreLossesOutput MagneticSimulator::calculate_core_losses(OperatingPoint& operatingPoint, Magnetic magnetic) {
    OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];
    if (!excitation.get_current()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Missing current in operating point");
    }

    double temperature = operatingPoint.get_conditions().get_ambient_temperature();
    double temperatureAfterLosses = temperature;
    SignalDescriptor magneticFluxDensity;
    CoreLossesOutput coreLossesOutput;

    do {
        temperature = temperatureAfterLosses;

        excitation = operatingPoint.get_excitations_per_winding()[0];
        operatingPoint.get_mutable_conditions().set_ambient_temperature(temperature);

        magneticFluxDensity = _magnetizingInductanceModel.calculate_inductance_and_magnetic_flux_density(magnetic.get_core(), magnetic.get_coil(), &operatingPoint).second;
        excitation.set_magnetic_flux_density(magneticFluxDensity);
        operatingPoint.get_mutable_excitations_per_winding()[0] = excitation;

        coreLossesOutput = _coreLossesModel.calculate_core_losses(magnetic.get_core(), excitation, temperature);

        // Core temperature from the full Temperature thermal-network model in core-only mode
        // (core nodes + convection/radiation, no turn/bobbin nodes) — more principled than the
        // legacy closed-form CoreTemperatureModel and still cheap (~0.07 ms/call).
        TemperatureConfig temperatureConfig;
        temperatureConfig.ambientTemperature = operatingPoint.get_conditions().get_ambient_temperature();
        temperatureConfig.coreLosses = coreLossesOutput.get_core_losses();
        temperatureConfig.coreOnly = true;
        temperatureConfig.plotSchematic = false;
        temperatureConfig.masCooling = operatingPoint.get_conditions().get_cooling();
        temperatureAfterLosses = Temperature(magnetic, temperatureConfig).calculateTemperatures().maximumTemperature;
        coreLossesOutput.set_temperature(temperatureAfterLosses);
    } while (fabs((temperature - temperatureAfterLosses) / temperatureAfterLosses) >= 0.05 && _enableTemperatureConvergence);


    return coreLossesOutput;
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

        if (inductanceOutput.get_leakage_inductance() &&
            !inductanceOutput.get_leakage_inductance()->get_leakage_inductance_per_winding().empty()) {
            electrical.set_leakage_inductance(
                inductanceOutput.get_leakage_inductance()->get_leakage_inductance_per_winding()[0]);
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
        double logStep = std::log10(stopFrequency / startFrequency) / (numberPoints - 1);

        std::vector<DatasheetImpedancePoint> impedancePoints;
        double maximumImpedanceMagnitude = std::numeric_limits<double>::lowest();
        for (size_t pointIndex = 0; pointIndex < numberPoints; ++pointIndex) {
            double frequency = startFrequency * std::pow(10.0, logStep * pointIndex);
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
        double coreLosses = outputs[operatingPointIndex].get_core_losses()->get_core_losses();
        double windingLosses = outputs[operatingPointIndex].get_winding_losses()
            ? outputs[operatingPointIndex].get_winding_losses()->get_winding_losses() : 0;

        TemperatureConfig config;
        config.ambientTemperature = ambient;
        config.masCooling = operatingPoints[operatingPointIndex].get_conditions().get_cooling();
        config.coreLosses = coreLosses;
        config.windingLosses = windingLosses;
        if (outputs[operatingPointIndex].get_winding_losses()) {
            config.windingLossesOutput = outputs[operatingPointIndex].get_winding_losses();
        }
        config.plotSchematic = false;

        double hotspot = Temperature(thermalMagnetic, config).calculateTemperatures().maximumTemperature;
        double temperatureRise = hotspot - ambient;
        if (temperatureRise > worstTemperatureRise) {
            worstTemperatureRise = temperatureRise;
            double totalLosses = coreLosses + windingLosses;
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
