#include "MagneticSimulator.h"
#include <MAS.hpp>


namespace OpenMagnetics {

MasWrapper MagneticSimulator::simulate(MasWrapper mas){
    return simulate(mas.get_mutable_inputs(), mas.get_mutable_magnetic());
}
MasWrapper MagneticSimulator::simulate(const InputsWrapper& inputs, const MagneticWrapper& magnetic){
    MasWrapper mas;
    std::vector<OutputsWrapper> outputs;
    std::vector<OperatingPoint> simulatedOperatingPoints;
    mas.set_inputs(inputs);

    for (auto& operatingPoint : mas.get_mutable_inputs().get_mutable_operating_points()){
        OutputsWrapper output;
        output.set_magnetizing_inductance(calculate_magnetizing_inductance(operatingPoint, magnetic));
        output.set_core_losses(calculate_core_loses(operatingPoint, magnetic));
        output.set_winding_losses(calculate_winding_losses(operatingPoint, magnetic, output.get_core_losses().value().get_temperature()));
        outputs.push_back(output);
        simulatedOperatingPoints.push_back(operatingPoint);
    }
    mas.get_mutable_inputs().set_operating_points(simulatedOperatingPoints);
    mas.set_magnetic(magnetic);
    mas.set_outputs(outputs);
    return mas;
}
MagnetizingInductanceOutput MagneticSimulator::calculate_magnetizing_inductance(OperatingPoint& operatingPoint, MagneticWrapper magnetic){
    return _magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil(), &operatingPoint);
}
WindingLossesOutput MagneticSimulator::calculate_winding_losses(OperatingPoint& operatingPoint, MagneticWrapper magnetic, std::optional<double> temperature){
    double simulationTemperature;
    if (!temperature) {
        simulationTemperature = operatingPoint.get_conditions().get_ambient_temperature();
    }
    else {
        simulationTemperature = temperature.value();
    }
    WindingLosses windingLosses;
    windingLosses.set_mirroring_dimension(0);
    windingLosses.set_winding_losses_harmonic_amplitude_threshold(0.05);
    return windingLosses.calculate_losses(magnetic, operatingPoint, simulationTemperature);
}

CoreLossesOutput MagneticSimulator::calculate_core_loses(OperatingPoint& operatingPoint, MagneticWrapper magnetic) {
    OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];
    if (!excitation.get_current()) {
        throw std::runtime_error("Missing current in operating point");
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

        coreLossesOutput = _coreLossesModel->get_core_losses(magnetic.get_core(), excitation, temperature);

        auto temperatureResult = _coreTemperatureModel->get_core_temperature(magnetic.get_core(), coreLossesOutput.get_core_losses(), operatingPoint.get_conditions().get_ambient_temperature());
        temperatureAfterLosses = temperatureResult.get_maximum_temperature();
        coreLossesOutput.set_temperature(temperatureAfterLosses);
    } while (fabs(temperature - temperatureAfterLosses) / temperatureAfterLosses >= 0.05 && _enableTemperatureConvergence);


    return coreLossesOutput;
}

} // namespace OpenMagnetics