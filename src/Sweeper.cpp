#include "Sweeper.h"
#include "ComplexPermeability.h"
#include "Constants.h"
#include "WindingLosses.h"
#include "CoreLosses.h"
#include "Settings.h"
#include <cmath>


namespace OpenMagnetics {


Curve2D Sweeper::sweep_impedance_over_frequency(MagneticWrapper magnetic, double start, double stop, size_t numberElements, std::string title) {
    // auto frequencies = linear_spaced_array(start, stop, numberElements);
    auto frequencies = logarithmic_spaced_array(start, stop, numberElements);

    std::vector<double> impedances;
    for (auto frequency : frequencies) {
        auto impedance = abs(OpenMagnetics::Impedance().calculate_impedance(magnetic, frequency));
        impedances.push_back(impedance);
    }

    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();
    return Curve2D(frequencies, impedances, title);
}

Curve2D Sweeper::sweep_resistance_over_frequency(MagneticWrapper magnetic, double start, double stop, size_t numberElements, size_t windingIndex, double temperature, std::string title) {
    auto frequencies = linear_spaced_array(start, stop, numberElements);
    // auto frequencies = logarithmic_spaced_array(start, stop, numberElements);

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

Curve2D Sweeper::sweep_core_resistance_over_frequency(MagneticWrapper magnetic, double start, double stop, size_t numberElements, double temperature, std::string title) {
    auto frequencies = linear_spaced_array(start, stop, numberElements);
    // auto frequencies = logarithmic_spaced_array(start, stop, numberElements);
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();

    auto magnetizingInductanceModel = MagnetizingInductance();
    auto magnetizingInductance = resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(core, coil).get_magnetizing_inductance());

    std::vector<double> coreResistances;
    auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(CoreLossesModels::STEINMETZ);
    for (auto frequency : frequencies) {
        auto coreResistance =  coreLossesModel->get_core_losses_series_resistance(core, frequency, temperature, magnetizingInductance);
        coreResistances.push_back(coreResistance);
        // std::cout << "frequency: " << frequency << std::endl;
        // std::cout << "coreResistance: " << coreResistance << std::endl;
    }

    return Curve2D(frequencies, coreResistances, title);
}

} // namespace OpenMagnetics
