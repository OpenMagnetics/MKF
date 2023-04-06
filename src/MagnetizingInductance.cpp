#include "MagnetizingInductance.h"

#include "../tests/TestingUtils.h"
#include "InputsWrapper.h"
#include "MagneticField.h"
#include "Reluctance.h"
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

void set_current_as_magnetizing_current(OperationPoint* operationPoint) {
    OperationPointExcitation excitation = InputsWrapper::get_primary_excitation(*operationPoint);
    auto current = excitation.get_current().value();

    auto currentExcitation = excitation.get_current().value();
    auto currentExcitationWaveform = currentExcitation.get_waveform().value();
    auto sampledCurrentWaveform =
        InputsWrapper::get_sampled_waveform(currentExcitationWaveform, excitation.get_frequency());

    if (sampledCurrentWaveform.get_data().size() > 0 && ((sampledCurrentWaveform.get_data().size() & (sampledCurrentWaveform.get_data().size() - 1)) != 0)) {
        throw std::invalid_argument("sampledCurrentWaveform vector size is not a power of 2");
    }

    currentExcitation.set_harmonics(
        InputsWrapper::get_harmonics_data(sampledCurrentWaveform, excitation.get_frequency()));
    currentExcitation.set_processed(InputsWrapper::get_processed_data(currentExcitation, sampledCurrentWaveform, true));
    excitation.set_current(currentExcitation);
    excitation.set_magnetizing_current(excitation.get_current().value());
    operationPoint->get_mutable_excitations_per_winding()[0] = excitation;
}

std::pair<double, ElectromagneticParameter> MagnetizingInductance::get_inductance_and_magnetic_flux_density(
    CoreWrapper core,
    WindingWrapper winding,
    OperationPoint* operationPoint) {

    InputsWrapper::make_waveform_size_power_of_two(operationPoint);

    std::pair<double, ElectromagneticParameter> result;
    double numberTurnsPrimary = winding.get_functional_description()[0].get_number_turns();
    double temperature = operationPoint->get_conditions().get_ambient_temperature();
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    double frequency = operationPoint->get_mutable_excitations_per_winding()[0].get_frequency();
    OperationPointExcitation excitation = InputsWrapper::get_primary_excitation(*operationPoint);
    OpenMagnetics::InitialPermeability initial_permeability;
    double currentInitialPermeability;
    double modifiedInitialPermeability;
    size_t timeout = 10;

    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(
        magic_enum::enum_cast<OpenMagnetics::ReluctanceModels>(_models["gapReluctance"]).value());
    double totalReluctance;
    double modifiedMagnetizingInductance = 5e-3;
    double currentMagnetizingInductance;

    if (!excitation.get_voltage()) {
        auto current = operationPoint->get_mutable_excitations_per_winding()[0].get_current().value();
        auto currentWaveform = current.get_waveform().value();
        if (!is_size_power_of_2(currentWaveform.get_data())) {
            auto currentSampledWaveform = InputsWrapper::get_sampled_waveform(currentWaveform, frequency);
            current.set_waveform(currentSampledWaveform);
            operationPoint->get_mutable_excitations_per_winding()[0].set_current(current);
        }
    }

    if (excitation.get_voltage()) {
        auto aux = operationPoint->get_mutable_excitations_per_winding()[0].get_voltage().value().get_waveform().value();
        if (aux.get_data().size() > 0 && ((aux.get_data().size() & (aux.get_data().size() - 1)) != 0)) {
            throw std::invalid_argument("voltage_data vector size is not a power of 2");
        }
    }
    if (excitation.get_current()) {
        auto aux = operationPoint->get_mutable_excitations_per_winding()[0].get_current().value().get_waveform().value();
        if (aux.get_data().size() > 0 && ((aux.get_data().size() & (aux.get_data().size() - 1)) != 0)) {
            throw std::invalid_argument("current_data vector size is not a power of 2");
        }
    }

    modifiedInitialPermeability = initial_permeability.get_initial_permeability(
        core.get_functional_description().get_material(), &temperature, nullptr, &frequency);
    if (!excitation.get_voltage()) {
        set_current_as_magnetizing_current(operationPoint);
        auto aux = operationPoint->get_mutable_excitations_per_winding()[0].get_magnetizing_current().value().get_waveform().value();
        if (aux.get_data().size() > 0 && ((aux.get_data().size() & (aux.get_data().size() - 1)) != 0)) {
            throw std::invalid_argument("magnetizing_current_data vector size is not a power of 2");
        }

    }

    do {
        currentMagnetizingInductance = modifiedMagnetizingInductance;


        do {
            currentInitialPermeability = modifiedInitialPermeability;
            totalReluctance = reluctanceModel->get_core_reluctance(core, currentInitialPermeability);
            modifiedMagnetizingInductance = pow(numberTurnsPrimary, 2) / totalReluctance;

            if (excitation.get_voltage()) {
                auto voltage = operationPoint->get_mutable_excitations_per_winding()[0].get_voltage().value();
                auto sampledVoltageWaveform =
                    InputsWrapper::get_sampled_waveform(voltage.get_waveform().value(), frequency);

                auto magnetizingCurrent = InputsWrapper::get_magnetizing_current(excitation, sampledVoltageWaveform,
                                                                                 modifiedMagnetizingInductance);

                auto sampledMagnetizingCurrentWaveform = InputsWrapper::get_sampled_waveform(
                    magnetizingCurrent.get_waveform().value(), excitation.get_frequency());
                magnetizingCurrent.set_harmonics(
                    InputsWrapper::get_harmonics_data(sampledMagnetizingCurrentWaveform, excitation.get_frequency()));
                magnetizingCurrent.set_processed(
                    InputsWrapper::get_processed_data(magnetizingCurrent, sampledMagnetizingCurrentWaveform, true));

                excitation.set_magnetizing_current(magnetizingCurrent);
                operationPoint->get_mutable_excitations_per_winding()[0] = excitation;
            }

            auto aux = operationPoint->get_mutable_excitations_per_winding()[0].get_magnetizing_current().value().get_waveform().value();
            if (aux.get_data().size() > 0 && ((aux.get_data().size() & (aux.get_data().size() - 1)) != 0)) {
                throw std::invalid_argument("magnetizing_current_data vector size is not a power of 2");
            }

            auto magneticFlux = OpenMagnetics::MagneticField::get_magnetic_flux(
                operationPoint->get_mutable_excitations_per_winding()[0].get_magnetizing_current().value(),
                totalReluctance, numberTurnsPrimary, frequency);
            auto magneticFluxDensity =
                OpenMagnetics::MagneticField::get_magnetic_flux_density(magneticFlux, effectiveArea, frequency);
            result.second = magneticFluxDensity;
            auto magneticFieldStrength = OpenMagnetics::MagneticField::get_magnetic_field_strength(
                magneticFluxDensity, currentInitialPermeability, frequency);

            modifiedInitialPermeability = initial_permeability.get_initial_permeability(
                core.get_functional_description().get_material(), &temperature,
                &magneticFieldStrength.get_processed().value().get_mutable_offset(), &frequency);

            timeout--;
            if (timeout == 0) {
                break;
            }
        } while (fabs(currentInitialPermeability - modifiedInitialPermeability) >= 1);

    } while (fabs(currentMagnetizingInductance - modifiedMagnetizingInductance) / modifiedMagnetizingInductance >=
             0.01);

    if (!excitation.get_voltage()) {
        operationPoint->get_mutable_excitations_per_winding()[0].set_voltage(
            InputsWrapper::get_induced_voltage(excitation, currentMagnetizingInductance));
    }

    result.first = currentMagnetizingInductance;
    return result;
}

double MagnetizingInductance::get_inductance_from_number_turns_and_gapping(CoreWrapper core,
                                                                           WindingWrapper winding,
                                                                           OperationPoint* operationPoint) {
    auto inductance_and_magnetic_flux_density = get_inductance_and_magnetic_flux_density(core, winding, operationPoint);

    return inductance_and_magnetic_flux_density.first;
}

int MagnetizingInductance::get_number_turns_from_gapping_and_inductance(CoreWrapper core, InputsWrapper* inputs) {
    auto operationPoint = inputs->get_operation_point(0);
    double desiredMagnetizingInductance =
        InputsWrapper::get_requirement_value(inputs->get_design_requirements().get_magnetizing_inductance());
    double temperature = operationPoint.get_conditions().get_ambient_temperature();
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    double frequency = operationPoint.get_mutable_excitations_per_winding()[0].get_frequency();
    OperationPointExcitation excitation = InputsWrapper::get_primary_excitation(operationPoint);
    OpenMagnetics::InitialPermeability initial_permeability;
    int numberTurnsPrimary;
    size_t timeout = 10;
    double currentInitialPermeability;
    double modifiedInitialPermeability;

    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(
        magic_enum::enum_cast<OpenMagnetics::ReluctanceModels>(_models["gapReluctance"]).value());
    double totalReluctance;

    currentInitialPermeability = initial_permeability.get_initial_permeability(
        core.get_functional_description().get_material(), &temperature, nullptr, &frequency);
    if (!excitation.get_voltage()) {
        set_current_as_magnetizing_current(&operationPoint);
        inputs->set_operation_point_by_index(operationPoint, 0);
    }
    while (true) {
        totalReluctance = reluctanceModel->get_core_reluctance(core, currentInitialPermeability);
        numberTurnsPrimary = std::round(sqrt(desiredMagnetizingInductance * totalReluctance));

        auto magneticFlux = OpenMagnetics::MagneticField::get_magnetic_flux(
            operationPoint.get_mutable_excitations_per_winding()[0].get_magnetizing_current().value(), totalReluctance,
            numberTurnsPrimary, frequency);
        auto magneticFluxDensity =
            OpenMagnetics::MagneticField::get_magnetic_flux_density(magneticFlux, effectiveArea, frequency);
        auto magneticFieldStrength = OpenMagnetics::MagneticField::get_magnetic_field_strength(
            magneticFluxDensity, currentInitialPermeability, frequency);

        modifiedInitialPermeability = initial_permeability.get_initial_permeability(
            core.get_functional_description().get_material(), &temperature,
            &magneticFieldStrength.get_processed().value().get_mutable_offset(), &frequency);

        if (fabs(currentInitialPermeability - modifiedInitialPermeability) < 1 || timeout == 0) {
            break;
        }
        else {
            currentInitialPermeability = modifiedInitialPermeability;
            timeout--;
        }
    }

    if (!excitation.get_voltage()) {
        operationPoint.get_mutable_excitations_per_winding()[0].set_voltage(
            InputsWrapper::get_induced_voltage(excitation, desiredMagnetizingInductance));
        inputs->set_operation_point_by_index(operationPoint, 0);
    }

    return std::max(1, numberTurnsPrimary);
}

CoreWrapper get_core_with_grinded_gapping(CoreWrapper core, double gapLength) {
    auto constants = OpenMagnetics::Constants();
    auto gappingJson = json::array();
    auto basicCentralGap = json();
    basicCentralGap["type"] = "subtractive";
    basicCentralGap["length"] = gapLength;
    auto basicLateralGap = json();
    basicLateralGap["type"] = "residual";
    basicLateralGap["length"] = constants.residualGap;
    std::vector<CoreGap> gapping;
    gapping.push_back(CoreGap(basicCentralGap));
    for (size_t i = 0; i < core.get_processed_description().value().get_columns().size() - 1; ++i) {
        gapping.push_back(CoreGap(basicLateralGap));
    }
    core.get_mutable_functional_description().set_gapping(gapping);
    core.process_gap();
    return core;
}

CoreWrapper get_core_with_distributed_gapping(CoreWrapper core, double gapLength, size_t numberDistributedGaps) {
    auto constants = OpenMagnetics::Constants();
    auto gappingJson = json::array();
    auto basicCentralGap = json();
    basicCentralGap["type"] = "subtractive";
    basicCentralGap["length"] = gapLength;
    auto basicLateralGap = json();
    basicLateralGap["type"] = "residual";
    basicLateralGap["length"] = constants.residualGap;
    std::vector<CoreGap> gapping;
    for (size_t i = 0; i < numberDistributedGaps; ++i) {
        gapping.push_back(CoreGap(basicCentralGap));
    }
    for (size_t i = 0; i < core.get_processed_description().value().get_columns().size() - 1; ++i) {
        gapping.push_back(CoreGap(basicLateralGap));
    }
    json coreJson;
    core.get_mutable_functional_description().set_gapping(gapping);
    core.process_gap();
    return core;
}

CoreWrapper get_core_with_spacer_gapping(CoreWrapper core, double gapLength) {
    auto gappingJson = json::array();
    auto basicCentralGap = json();
    basicCentralGap["type"] = "additive";
    basicCentralGap["length"] = gapLength;
    auto basicLateralGap = json();
    basicLateralGap["type"] = "additive";
    basicLateralGap["length"] = gapLength;
    std::vector<CoreGap> gapping;
    gapping.push_back(CoreGap(basicCentralGap));
    for (size_t i = 0; i < core.get_processed_description().value().get_columns().size() - 1; ++i) {
        gapping.push_back(CoreGap(basicLateralGap));
    }
    core.get_mutable_functional_description().set_gapping(gapping);
    core.process_gap();
    return core;
}

std::vector<CoreGap> MagnetizingInductance::get_gapping_from_number_turns_and_inductance(CoreWrapper core,
                                                                                         WindingWrapper winding,
                                                                                         InputsWrapper* inputs,
                                                                                         GappingType gappingType,
                                                                                         size_t decimals) {
    auto constants = OpenMagnetics::Constants();
    auto operationPoint = inputs->get_operation_point(0);
    double numberTurnsPrimary = winding.get_functional_description()[0].get_number_turns();
    double desiredMagnetizingInductance =
        InputsWrapper::get_requirement_value(inputs->get_design_requirements().get_magnetizing_inductance());
    double temperature = operationPoint.get_conditions().get_ambient_temperature();
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    OperationPointExcitation excitation = InputsWrapper::get_primary_excitation(operationPoint);
    double frequency = operationPoint.get_mutable_excitations_per_winding()[0].get_frequency();
    OpenMagnetics::InitialPermeability initial_permeability;
    size_t timeout = 10;
    double modifiedInitialPermeability;
    double currentInitialPermeability;

    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(
        magic_enum::enum_cast<OpenMagnetics::ReluctanceModels>(_models["gapReluctance"]).value());
    double neededTotalReluctance = pow(numberTurnsPrimary, 2) / desiredMagnetizingInductance;

    currentInitialPermeability = initial_permeability.get_initial_permeability(
        core.get_functional_description().get_material(), &temperature, nullptr, &frequency);
    if (!excitation.get_voltage()) {
        set_current_as_magnetizing_current(&operationPoint);
        inputs->set_operation_point_by_index(operationPoint, 0);
    }

    while (true) {
        auto magneticFlux = OpenMagnetics::MagneticField::get_magnetic_flux(
            operationPoint.get_mutable_excitations_per_winding()[0].get_magnetizing_current().value(),
            neededTotalReluctance, numberTurnsPrimary, frequency);
        auto magneticFluxDensity =
            OpenMagnetics::MagneticField::get_magnetic_flux_density(magneticFlux, effectiveArea, frequency);
        auto magneticFieldStrength = OpenMagnetics::MagneticField::get_magnetic_field_strength(
            magneticFluxDensity, currentInitialPermeability, frequency);

        modifiedInitialPermeability = initial_permeability.get_initial_permeability(
            core.get_functional_description().get_material(), &temperature,
            &magneticFieldStrength.get_processed().value().get_mutable_offset(), &frequency);

        if (fabs(currentInitialPermeability - modifiedInitialPermeability) < 1 || timeout == 0) {
            break;
        }
        else {
            currentInitialPermeability = modifiedInitialPermeability;
            timeout--;
        }
    }

    if (!excitation.get_voltage()) {
        operationPoint.get_mutable_excitations_per_winding()[0].set_voltage(
            InputsWrapper::get_induced_voltage(excitation, desiredMagnetizingInductance));
        inputs->set_operation_point_by_index(operationPoint, 0);
    }

    double gapLength = constants.residualGap;
    double gapLengthModification = constants.initial_gap_length_for_searching;
    bool increasingGap = true;
    double fringingFactorOneGap = 0;

    double reluctance = 0;
    size_t numberDistributedGaps = 3;
    timeout = 100;
    CoreWrapper gappedCore;

    while (true) {
        reluctance = 0;
        switch (gappingType) {
            case GappingType::GRINDED:
                gappedCore = get_core_with_grinded_gapping(core, gapLength);
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
                    fringingFactorOneGap =
                        reluctanceModel->get_gap_reluctance(gappedCore.get_gapping()[0])["fringing_factor"];
                    if (fringingFactorOneGap < constants.minimum_distributed_fringing_factor &&
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
                    fringingFactorOneGap =
                        reluctanceModel->get_gap_reluctance(gappedCore.get_gapping()[0])["fringing_factor"];
                    if (fringingFactorOneGap > constants.maximum_distributed_fringing_factor) {
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
                throw std::runtime_error("Unknown type of gap, options are {GRINDED, SPACER, RESIDUAL, DISTRIBUTED}");
        }

        reluctance = reluctanceModel->get_core_reluctance(gappedCore, currentInitialPermeability);

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
        case GappingType::GRINDED:
            return get_core_with_grinded_gapping(core, gapLength).get_gapping();
        case GappingType::SPACER:
            return get_core_with_spacer_gapping(core, gapLength).get_gapping();
        case GappingType::RESIDUAL:
            throw std::runtime_error("Residual type cannot be chosen to calculate the needed gapping");
            break;
        case GappingType::DISTRIBUTED:
            return get_core_with_distributed_gapping(core, gapLength, numberDistributedGaps).get_gapping();
            break;
        default:
            throw std::runtime_error("Unknown type of gap, options are {GRINDED, SPACER, RESIDUAL, DISTRIBUTED}");
    }
}

} // namespace OpenMagnetics
