#include "MagnetizingInductance.h"

#include "InputsWrapper.h"
#include "MagneticField.h"
#include "Reluctance.h"
#include "Settings.h"
#include "Utils.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>

namespace OpenMagnetics {

double calculate_air_inductance(int64_t numberTurnsPrimary, CoreWrapper core) {

    auto bobbin = BobbinWrapper::create_quick_bobbin(core);
    auto bobbinColumnDepth = bobbin.get_processed_description()->get_column_depth();
    auto bobbinColumnWidth = bobbin.get_processed_description()->get_column_width().value();
    auto bobbinWindingWindowWidth = bobbin.get_winding_window_dimensions()[0];
    auto bobbinWindingWindowHeight = bobbin.get_winding_window_dimensions()[1];
    auto meanLengthRadius = (bobbinColumnDepth + bobbinColumnWidth) / 2 + bobbinWindingWindowWidth / 4;

    double coilInternalArea = std::numbers::pi * pow(meanLengthRadius, 2);

    double coreColumnArea = core.get_processed_description()->get_columns()[0].get_area();

    double airAreaProportion = (coilInternalArea - coreColumnArea) / coilInternalArea;

    return Constants().vacuumPermeability * pow(numberTurnsPrimary, 2) * (coilInternalArea * airAreaProportion * 2) / bobbinWindingWindowHeight;
}

std::pair<MagnetizingInductanceOutput, SignalDescriptor> MagnetizingInductance::calculate_inductance_and_magnetic_flux_density(
    CoreWrapper core,
    CoilWrapper winding,
    OperatingPoint* operatingPoint) {

    double frequency = Defaults().coreAdviserFrequencyReference;
    double temperature = Defaults().ambientTemperature;

    if (operatingPoint) {
        InputsWrapper::make_waveform_size_power_of_two(operatingPoint);
        temperature = operatingPoint->get_conditions().get_ambient_temperature();
        frequency = operatingPoint->get_mutable_excitations_per_winding()[0].get_frequency();
    }

    std::pair<MagnetizingInductanceOutput, SignalDescriptor> result;
    double numberWindings = winding.get_functional_description().size();
    double numberTurnsPrimary = winding.get_functional_description()[0].get_number_turns();
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    OpenMagnetics::InitialPermeability initialPermeability;
    double currentInitialPermeability;
    double modifiedInitialPermeability;
    size_t timeout = 10;

    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(magic_enum::enum_cast<OpenMagnetics::ReluctanceModels>(_models["gapReluctance"]).value());
    double totalReluctance;
    double modifiedMagnetizingInductance = 5e-3;
    double currentMagnetizingInductance;

    if (operatingPoint) {
        OperatingPointExcitation excitation = InputsWrapper::get_primary_excitation(*operatingPoint);
        if (!excitation.get_voltage()) {
            auto current = operatingPoint->get_mutable_excitations_per_winding()[0].get_current().value();
            auto currentWaveform = current.get_waveform().value();
            if (!is_size_power_of_2(currentWaveform.get_data())) {
                auto currentSampledWaveform = InputsWrapper::calculate_sampled_waveform(currentWaveform, frequency);
                current.set_waveform(currentSampledWaveform);
                operatingPoint->get_mutable_excitations_per_winding()[0].set_current(current);
            }
        }

        if (excitation.get_voltage()) {
            auto aux = operatingPoint->get_mutable_excitations_per_winding()[0].get_voltage().value().get_waveform().value();
            if (aux.get_data().size() > 0 && ((aux.get_data().size() & (aux.get_data().size() - 1)) != 0)) {
                throw std::invalid_argument("voltage_data vector size is not a power of 2");
            }
        }
        if (excitation.get_current()) {
            auto aux = operatingPoint->get_mutable_excitations_per_winding()[0].get_current().value().get_waveform().value();
            if (aux.get_data().size() > 0 && ((aux.get_data().size() & (aux.get_data().size() - 1)) != 0)) {
                throw std::invalid_argument("current_data vector size is not a power of 2");
            }
        }
        if (!excitation.get_voltage()) {
            InputsWrapper::set_current_as_magnetizing_current(operatingPoint);
            auto aux = operatingPoint->get_mutable_excitations_per_winding()[0].get_magnetizing_current().value().get_waveform().value();
            if (aux.get_data().size() > 0 && ((aux.get_data().size() & (aux.get_data().size() - 1)) != 0)) {
                throw std::invalid_argument("magnetizing_current_data vector size from current is not a power of 2");
            }

        }
    }

    modifiedInitialPermeability = initialPermeability.get_initial_permeability(core.resolve_material(), temperature, std::nullopt, frequency);

    MagnetizingInductanceOutput magnetizingInductanceOutput;

    do {
        currentMagnetizingInductance = modifiedMagnetizingInductance;


        do {
            currentInitialPermeability = modifiedInitialPermeability;
            magnetizingInductanceOutput = reluctanceModel->get_core_reluctance(core, currentInitialPermeability);
            totalReluctance = magnetizingInductanceOutput.get_core_reluctance();
            modifiedMagnetizingInductance = pow(numberTurnsPrimary, 2) / totalReluctance;


            if (operatingPoint) {
                OperatingPointExcitation excitation = InputsWrapper::get_primary_excitation(*operatingPoint);
                if (numberWindings == 1 && excitation.get_current()) {
                    InputsWrapper::set_current_as_magnetizing_current(operatingPoint);
                }
                else if (InputsWrapper::is_multiport_inductor(*operatingPoint)) {
                    auto magnetizingCurrent = InputsWrapper::get_multiport_inductor_magnetizing_current(*operatingPoint);
                    excitation.set_magnetizing_current(magnetizingCurrent);
                    operatingPoint->get_mutable_excitations_per_winding()[0] = excitation;
                }
                else if (excitation.get_voltage()) {
                    auto voltage = operatingPoint->get_mutable_excitations_per_winding()[0].get_voltage().value();
                    auto sampledVoltageWaveform = InputsWrapper::calculate_sampled_waveform(voltage.get_waveform().value(), frequency);

                    auto magnetizingCurrent = InputsWrapper::calculate_magnetizing_current(excitation,
                                                                                           sampledVoltageWaveform,
                                                                                           modifiedMagnetizingInductance,
                                                                                           false);

                    auto sampledMagnetizingCurrentWaveform = InputsWrapper::calculate_sampled_waveform(magnetizingCurrent.get_waveform().value(), excitation.get_frequency());
                    magnetizingCurrent.set_harmonics(InputsWrapper::calculate_harmonics_data(sampledMagnetizingCurrentWaveform, excitation.get_frequency()));
                    magnetizingCurrent.set_processed(InputsWrapper::calculate_processed_data(magnetizingCurrent, sampledMagnetizingCurrentWaveform, false));

                    excitation.set_magnetizing_current(magnetizingCurrent);
                    operatingPoint->get_mutable_excitations_per_winding()[0] = excitation;
                }

                auto aux = operatingPoint->get_mutable_excitations_per_winding()[0].get_magnetizing_current().value().get_waveform().value();
                if (aux.get_data().size() > 0 && ((aux.get_data().size() & (aux.get_data().size() - 1)) != 0)) {
                    throw std::invalid_argument("magnetizing_current_data vector size from voltage is not a power of 2");
                }

                auto magneticFlux = OpenMagnetics::MagneticField::calculate_magnetic_flux(operatingPoint->get_mutable_excitations_per_winding()[0].get_magnetizing_current().value(), totalReluctance, numberTurnsPrimary);
                auto magneticFluxDensity = OpenMagnetics::MagneticField::calculate_magnetic_flux_density(magneticFlux, effectiveArea);
                result.second = magneticFluxDensity;
                auto magneticFieldStrength = OpenMagnetics::MagneticField::calculate_magnetic_field_strength(magneticFluxDensity, currentInitialPermeability);

                modifiedInitialPermeability = initialPermeability.get_initial_permeability(core.resolve_material(), temperature, magneticFieldStrength.get_processed().value().get_mutable_offset(), frequency);
            }

            timeout--;
            if (timeout == 0) {
                break;
            }
        } while (fabs(currentInitialPermeability - modifiedInitialPermeability) >= 1);

    } while (fabs(currentMagnetizingInductance - modifiedMagnetizingInductance) / modifiedMagnetizingInductance >= 0.01);

    if (operatingPoint) {
        OperatingPointExcitation excitation = InputsWrapper::get_primary_excitation(*operatingPoint);
        if (!excitation.get_voltage()) {
            operatingPoint->get_mutable_excitations_per_winding()[0].set_voltage(InputsWrapper::calculate_induced_voltage(excitation, currentMagnetizingInductance));
        }
    }

    auto settings = Settings::GetInstance();
    if (settings->get_magnetizing_inductance_include_air_inductance()) {
        currentMagnetizingInductance += calculate_air_inductance(numberTurnsPrimary, core);
    }

    DimensionWithTolerance magnetizingInductanceWithTolerance;
    magnetizingInductanceWithTolerance.set_nominal(currentMagnetizingInductance);
    magnetizingInductanceOutput.set_magnetizing_inductance(magnetizingInductanceWithTolerance);

    result.first = magnetizingInductanceOutput;
    return result;
}

MagnetizingInductanceOutput MagnetizingInductance::calculate_inductance_from_number_turns_and_gapping(CoreWrapper core,
                                                                           CoilWrapper winding,
                                                                           OperatingPoint* operatingPoint) {
    auto inductance_and_magnetic_flux_density = calculate_inductance_and_magnetic_flux_density(core, winding, operatingPoint);

    return inductance_and_magnetic_flux_density.first;
}

int MagnetizingInductance::calculate_number_turns_from_gapping_and_inductance(CoreWrapper core, InputsWrapper* inputs, DimensionalValues preferredValue) {
    double desiredMagnetizingInductance = resolve_dimensional_values(inputs->get_design_requirements().get_magnetizing_inductance(), preferredValue);
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    double frequency = Defaults().coreAdviserFrequencyReference;
    double temperature = Defaults().ambientTemperature;
    OperatingPoint operatingPoint;

    if (inputs->get_operating_points().size() > 0) {
        operatingPoint = inputs->get_operating_point(0);
        temperature = operatingPoint.get_conditions().get_ambient_temperature();
        frequency = operatingPoint.get_mutable_excitations_per_winding()[0].get_frequency();
    }
    OpenMagnetics::InitialPermeability initialPermeability;
    int numberTurnsPrimary;
    size_t timeout = 10;
    double currentInitialPermeability;
    double modifiedInitialPermeability;

    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(magic_enum::enum_cast<OpenMagnetics::ReluctanceModels>(_models["gapReluctance"]).value());
    double totalReluctance;

    currentInitialPermeability = initialPermeability.get_initial_permeability(core.resolve_material(), temperature, std::nullopt, frequency);
    if (inputs->get_operating_points().size() > 0) {
        OperatingPointExcitation excitation = InputsWrapper::get_primary_excitation(operatingPoint);
        if (!excitation.get_voltage()) {
            InputsWrapper::set_current_as_magnetizing_current(&operatingPoint);
            inputs->set_operating_point_by_index(operatingPoint, 0);
        }
    }

    while (true) {

        auto magnetizingInductanceOutput = reluctanceModel->get_core_reluctance(core, currentInitialPermeability);
        totalReluctance = magnetizingInductanceOutput.get_core_reluctance();
        numberTurnsPrimary = std::round(sqrt(desiredMagnetizingInductance * totalReluctance));

        if (inputs->get_operating_points().size() > 0) {
            auto magneticFlux = OpenMagnetics::MagneticField::calculate_magnetic_flux(operatingPoint.get_mutable_excitations_per_winding()[0].get_magnetizing_current().value(), totalReluctance, numberTurnsPrimary);
            auto magneticFluxDensity = OpenMagnetics::MagneticField::calculate_magnetic_flux_density(magneticFlux, effectiveArea);
            auto magneticFieldStrength = OpenMagnetics::MagneticField::calculate_magnetic_field_strength(magneticFluxDensity, currentInitialPermeability);

            modifiedInitialPermeability = initialPermeability.get_initial_permeability(core.resolve_material(), temperature, magneticFieldStrength.get_processed().value().get_mutable_offset(), frequency);

            if (fabs(currentInitialPermeability - modifiedInitialPermeability) < 1 || timeout == 0) {
                break;
            }
            else {
                currentInitialPermeability = modifiedInitialPermeability;
                timeout--;
            }
        }
    }

    if (inputs->get_operating_points().size() > 0) {
        OperatingPointExcitation excitation = InputsWrapper::get_primary_excitation(operatingPoint);
        if (!excitation.get_voltage()) {
            operatingPoint.get_mutable_excitations_per_winding()[0].set_voltage(InputsWrapper::calculate_induced_voltage(excitation, desiredMagnetizingInductance));
            inputs->set_operating_point_by_index(operatingPoint, 0);
        }
    }

    return std::max(1, numberTurnsPrimary);
}

CoreWrapper get_core_with_ground_gapping(CoreWrapper core, double gapLength) {
    auto constants = OpenMagnetics::Constants();
    auto basicCentralGap = CoreGap();
    basicCentralGap.set_type(GapType::SUBTRACTIVE);
    basicCentralGap.set_length(gapLength);
    auto basicLateralGap = CoreGap();
    basicLateralGap.set_type(GapType::RESIDUAL);
    basicLateralGap.set_length(constants.residualGap);
    std::vector<CoreGap> gapping;
    gapping.push_back(basicCentralGap);
    for (size_t i = 0; i < core.get_processed_description().value().get_columns().size() - 1; ++i) {
        gapping.push_back(basicLateralGap);
    }
    core.get_mutable_functional_description().set_gapping(gapping);
    core.process_gap();
    return core;
}

CoreWrapper get_core_with_distributed_gapping(CoreWrapper core, double gapLength, size_t numberDistributedGaps) {
    auto constants = OpenMagnetics::Constants();
    auto basicCentralGap = CoreGap();
    basicCentralGap.set_type(GapType::SUBTRACTIVE);
    basicCentralGap.set_length(gapLength);
    auto basicLateralGap = CoreGap();
    basicLateralGap.set_type(GapType::RESIDUAL);
    basicLateralGap.set_length(constants.residualGap);
    std::vector<CoreGap> gapping;
    for (size_t i = 0; i < numberDistributedGaps; ++i) {
        gapping.push_back(basicCentralGap);
    }
    for (size_t i = 0; i < core.get_processed_description().value().get_columns().size() - 1; ++i) {
        gapping.push_back(basicLateralGap);
    }
    core.get_mutable_functional_description().set_gapping(gapping);
    core.process_gap();
    return core;
}

CoreWrapper get_core_with_spacer_gapping(CoreWrapper core, double gapLength) {
    auto basicCentralGap = CoreGap();
    basicCentralGap.set_type(GapType::ADDITIVE);
    basicCentralGap.set_length(gapLength);
    auto basicLateralGap = CoreGap();
    basicLateralGap.set_type(GapType::ADDITIVE);
    basicLateralGap.set_length(gapLength);
    std::vector<CoreGap> gapping;
    gapping.push_back(basicCentralGap);
    for (size_t i = 0; i < core.get_processed_description().value().get_columns().size() - 1; ++i) {
        gapping.push_back(basicLateralGap);
    }
    core.get_mutable_functional_description().set_gapping(gapping);
    core.process_gap();
    return core;
}

std::vector<CoreGap> MagnetizingInductance::calculate_gapping_from_number_turns_and_inductance(CoreWrapper core,
                                                                                         CoilWrapper winding,
                                                                                         InputsWrapper* inputs,
                                                                                         GappingType gappingType,
                                                                                         size_t decimals) {
    double frequency = Defaults().coreAdviserFrequencyReference;
    double temperature = Defaults().ambientTemperature;
    OperatingPointExcitation excitation;
    OperatingPoint operatingPoint;

    auto constants = OpenMagnetics::Constants();
    if (inputs->get_operating_points().size() > 0) {
        operatingPoint = inputs->get_operating_point(0);
        excitation = InputsWrapper::get_primary_excitation(operatingPoint);
        temperature = operatingPoint.get_conditions().get_ambient_temperature();
        frequency = operatingPoint.get_mutable_excitations_per_winding()[0].get_frequency();
    }

    double numberTurnsPrimary = winding.get_functional_description()[0].get_number_turns();
    double desiredMagnetizingInductance = resolve_dimensional_values(inputs->get_design_requirements().get_magnetizing_inductance());
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    OpenMagnetics::InitialPermeability initialPermeability;
    size_t timeout = 10;
    double modifiedInitialPermeability;
    double currentInitialPermeability;

    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(magic_enum::enum_cast<OpenMagnetics::ReluctanceModels>(_models["gapReluctance"]).value());
    double neededTotalReluctance = pow(numberTurnsPrimary, 2) / desiredMagnetizingInductance;

    currentInitialPermeability = initialPermeability.get_initial_permeability(core.resolve_material(), temperature, std::nullopt, frequency);

    if (!excitation.get_voltage() && excitation.get_current()) {
        InputsWrapper::set_current_as_magnetizing_current(&operatingPoint);
        inputs->set_operating_point_by_index(operatingPoint, 0);
    }

    while (true) {
        if (!excitation.get_magnetizing_current()) {
            break;
        }
        auto magneticFlux = OpenMagnetics::MagneticField::calculate_magnetic_flux(operatingPoint.get_mutable_excitations_per_winding()[0].get_magnetizing_current().value(), neededTotalReluctance, numberTurnsPrimary);
        auto magneticFluxDensity = OpenMagnetics::MagneticField::calculate_magnetic_flux_density(magneticFlux, effectiveArea);
        auto magneticFieldStrength = OpenMagnetics::MagneticField::calculate_magnetic_field_strength(magneticFluxDensity, currentInitialPermeability);

        modifiedInitialPermeability = initialPermeability.get_initial_permeability(core.resolve_material(), temperature, magneticFieldStrength.get_processed().value().get_mutable_offset(), frequency);

        if (fabs(currentInitialPermeability - modifiedInitialPermeability) < 1 || timeout == 0) {
            break;
        }
        else {
            currentInitialPermeability = modifiedInitialPermeability;
            timeout--;
        }
    }

    if (!excitation.get_voltage() && excitation.get_current()) {
        operatingPoint.get_mutable_excitations_per_winding()[0].set_voltage(InputsWrapper::calculate_induced_voltage(excitation, desiredMagnetizingInductance));
        inputs->set_operating_point_by_index(operatingPoint, 0);
    }

    double gapLength = constants.residualGap;
    double gapLengthModification = constants.initialGapLengthForSearching;
    bool increasingGap = true;
    double fringingFactorOneGap = 0;

    double reluctance = 0;
    size_t numberDistributedGaps = 3;
    timeout = 100;
    CoreWrapper gappedCore;

    while (true) {
        reluctance = 0;
        switch (gappingType) {
            case GappingType::GROUND:
                gappedCore = get_core_with_ground_gapping(core, gapLength);
                break;
            case GappingType::SPACER:
                gappedCore = get_core_with_spacer_gapping(core, gapLength);
                break;
            case GappingType::RESIDUAL:
                throw std::runtime_error("Residual type cannot be chosen to calculate the needed gapping");
                break;
            case GappingType::DISTRIBUTED:
                while (numberDistributedGaps > 3) {
                    gappedCore = get_core_with_distributed_gapping(core, gapLength, numberDistributedGaps);
                    fringingFactorOneGap = reluctanceModel->get_gap_reluctance(gappedCore.get_gapping()[0]).get_fringing_factor();
                    if (fringingFactorOneGap < constants.minimumDistributedFringingFactor &&
                        numberDistributedGaps > 1) {
                        gapLength *= numberDistributedGaps;
                        numberDistributedGaps -= 2;
                        gapLength /= numberDistributedGaps;
                    }
                    else {
                        break;
                    }
                }
                while (true) {
                    gappedCore = get_core_with_distributed_gapping(core, gapLength, numberDistributedGaps);
                    fringingFactorOneGap = reluctanceModel->get_gap_reluctance(gappedCore.get_gapping()[0]).get_fringing_factor();
                    if (fringingFactorOneGap > constants.maximumDistributedFringingFactor) {
                        gapLength *= numberDistributedGaps;
                        numberDistributedGaps += 2;
                        gapLength /= numberDistributedGaps;
                    }
                    else {
                        break;
                    }
                }
                break;
            default:
                throw std::runtime_error("Unknown type of gap, options are {GROUND, SPACER, RESIDUAL, DISTRIBUTED}");
        }

        auto magnetizingInductanceOutput = reluctanceModel->get_core_reluctance(gappedCore, currentInitialPermeability);
        reluctance = magnetizingInductanceOutput.get_core_reluctance();

        if (fabs(neededTotalReluctance - reluctance) / neededTotalReluctance < 0.001 || timeout == 0) {
            break;
        }
        else {
            if (neededTotalReluctance < reluctance && increasingGap) {
                increasingGap = false;
                gapLengthModification = std::max(gapLengthModification / 2., constants.residualGap);
            }
            if (neededTotalReluctance > reluctance && !increasingGap) {
                increasingGap = true;
                gapLengthModification = std::max(gapLengthModification / 2., constants.residualGap);
            }
            if (increasingGap) {
                gapLength += gapLengthModification;
            }
            else {
                gapLength -= gapLengthModification;
            }

            timeout--;
        }
    }

    gapLength = roundFloat(gapLength, decimals);

    switch (gappingType) {
        case GappingType::GROUND:
            return get_core_with_ground_gapping(core, gapLength).get_gapping();
        case GappingType::SPACER:
            return get_core_with_spacer_gapping(core, gapLength).get_gapping();
        case GappingType::RESIDUAL:
            throw std::runtime_error("Residual type cannot be chosen to calculate the needed gapping");
            break;
        case GappingType::DISTRIBUTED:
            return get_core_with_distributed_gapping(core, gapLength, numberDistributedGaps).get_gapping();
            break;
        default:
            throw std::runtime_error("Unknown type of gap, options are {GROUND, SPACER, RESIDUAL, DISTRIBUTED}");
    }
}

} // namespace OpenMagnetics
