#include <fstream>
#include <numbers>
#include <iostream>
#include <cmath>
#include <vector>
#include <filesystem>
#include <streambuf>
#include <magic_enum.hpp>
#include "MagnetizingInductance.h"
#include "MagneticField.h"
#include "Reluctance.h"
#include "Utils.h"
// #include "../tests/TestingUtils.h"


namespace OpenMagnetics {

    double MagnetizingInductance::get_inductance_from_number_turns_and_gapping(CoreWrapper core,
                                                                               WindingWrapper winding,
                                                                               OperationPoint operationPoint){

        double numberTurnsPrimary = winding.get_functional_description()[0].get_number_turns();
        double temperature = operationPoint.get_conditions().get_ambient_temperature();
        double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
        double frequency = operationPoint.get_mutable_excitations_per_winding()[0].get_frequency();
        OpenMagnetics::InitialPermeability initial_permeability;
        double currentInitialPermeability;
        double modifiedInitialPermeability;
        size_t timeout = 10;

        auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(magic_enum::enum_cast<OpenMagnetics::ReluctanceModels>(_models["gapReluctance"]).value());
        double totalReluctance;
        double modifiedMagnetizingInductance = 5e-3;
        double currentMagnetizingInductance;

        modifiedInitialPermeability = initial_permeability.get_initial_permeability(core.get_functional_description().get_material(),
                                                                                    &temperature,
                                                                                    nullptr,
                                                                                    &frequency);

        do {
            currentMagnetizingInductance = modifiedMagnetizingInductance;


            if (!operationPoint.get_mutable_excitations_per_winding()[0].get_magnetizing_current()) {
                operationPoint = InputsWrapper::process_operation_point(operationPoint, currentMagnetizingInductance);
            }

            do
            {
                currentInitialPermeability = modifiedInitialPermeability;
                totalReluctance = reluctanceModel->get_core_reluctance(core, currentInitialPermeability);

                auto magneticFlux = OpenMagnetics::MagneticField::get_magnetic_flux(operationPoint.get_mutable_excitations_per_winding()[0].get_magnetizing_current().value(),
                                                                                    totalReluctance,
                                                                                    numberTurnsPrimary,
                                                                                    frequency);
                auto magneticFluxDensity = OpenMagnetics::MagneticField::get_magnetic_flux_density(magneticFlux, effectiveArea, frequency);
                auto magneticFieldStrength = OpenMagnetics::MagneticField::get_magnetic_field_strength(magneticFluxDensity, currentInitialPermeability, frequency);

                modifiedInitialPermeability = initial_permeability.get_initial_permeability(core.get_functional_description().get_material(),
                                                                                    &temperature,
                                                                                    &magneticFieldStrength.get_processed().value().get_mutable_offset(),
                                                                                    &frequency);

                timeout--;
                if (timeout == 0) {
                    break;
                }
            } while (fabs(currentInitialPermeability - modifiedInitialPermeability) >= 1);

            modifiedMagnetizingInductance = pow(numberTurnsPrimary, 2) / totalReluctance;
        } while (fabs(currentMagnetizingInductance - modifiedMagnetizingInductance) / modifiedMagnetizingInductance >= 0.01);

        return currentMagnetizingInductance;
    }

    int MagnetizingInductance::get_number_turns_from_gapping_and_inductance(CoreWrapper core,
                                                                               InputsWrapper inputs){
        auto operationPoint = inputs.get_operation_point(0);
        double desiredMagnetizingInductance = InputsWrapper::get_requirement_value(inputs.get_design_requirements().get_magnetizing_inductance());
        double temperature = operationPoint.get_conditions().get_ambient_temperature();
        double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
        double frequency = operationPoint.get_mutable_excitations_per_winding()[0].get_frequency();
        OpenMagnetics::InitialPermeability initial_permeability;
        int numberTurnsPrimary;
        size_t timeout = 10;
        double currentInitialPermeability;
        double modifiedInitialPermeability;

        auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(magic_enum::enum_cast<OpenMagnetics::ReluctanceModels>(_models["gapReluctance"]).value());
        double totalReluctance;

        currentInitialPermeability = initial_permeability.get_initial_permeability(core.get_functional_description().get_material(),
                                                                                   &temperature,
                                                                                   nullptr,
                                                                                   &frequency);

        while (true) {

            totalReluctance = reluctanceModel->get_core_reluctance(core, currentInitialPermeability);
            numberTurnsPrimary =  std::round(sqrt(desiredMagnetizingInductance * totalReluctance));

            auto magneticFlux = OpenMagnetics::MagneticField::get_magnetic_flux(operationPoint.get_mutable_excitations_per_winding()[0].get_magnetizing_current().value(),
                                                                                totalReluctance,
                                                                                numberTurnsPrimary,
                                                                                frequency);
            auto magneticFluxDensity = OpenMagnetics::MagneticField::get_magnetic_flux_density(magneticFlux, effectiveArea, frequency);
            auto magneticFieldStrength = OpenMagnetics::MagneticField::get_magnetic_field_strength(magneticFluxDensity, currentInitialPermeability, frequency);

            modifiedInitialPermeability = initial_permeability.get_initial_permeability(core.get_functional_description().get_material(),
                                                                                        &temperature,
                                                                                        &magneticFieldStrength.get_processed().value().get_mutable_offset(),
                                                                                        &frequency);

            if (fabs(currentInitialPermeability - modifiedInitialPermeability) < 1 || timeout == 0) {
                break;
            }
            else {
                currentInitialPermeability = modifiedInitialPermeability;
                timeout--;
            }
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
        for (size_t i=0; i<core.get_processed_description().value().get_columns().size() - 1; ++i) {
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
        for (size_t i=0; i<numberDistributedGaps; ++i) {
            gapping.push_back(CoreGap(basicCentralGap));
        }
        for (size_t i=0; i<core.get_processed_description().value().get_columns().size() - 1; ++i) {
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
        for (size_t i=0; i<core.get_processed_description().value().get_columns().size() - 1; ++i) {
            gapping.push_back(CoreGap(basicLateralGap));
        }
        core.get_mutable_functional_description().set_gapping(gapping);
        core.process_gap();
        return core;
    }

    std::vector<CoreGap> MagnetizingInductance::get_gapping_from_number_turns_and_inductance(CoreWrapper core,
                                                                              WindingWrapper winding,
                                                                              InputsWrapper inputs,
                                                                              GappingType gappingType, 
                                                                              size_t decimals){
        auto constants = OpenMagnetics::Constants();
        auto operationPoint = inputs.get_operation_point(0);
        double numberTurnsPrimary = winding.get_functional_description()[0].get_number_turns();
        double desiredMagnetizingInductance = InputsWrapper::get_requirement_value(inputs.get_design_requirements().get_magnetizing_inductance());
        double temperature = operationPoint.get_conditions().get_ambient_temperature();
        double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
        double frequency = operationPoint.get_mutable_excitations_per_winding()[0].get_frequency();
        OpenMagnetics::InitialPermeability initial_permeability;
        size_t timeout = 10;
        double modifiedInitialPermeability;
        double currentInitialPermeability;

        auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(magic_enum::enum_cast<OpenMagnetics::ReluctanceModels>(_models["gapReluctance"]).value());
        double neededTotalReluctance = pow(numberTurnsPrimary, 2) / desiredMagnetizingInductance;

        currentInitialPermeability = initial_permeability.get_initial_permeability(core.get_functional_description().get_material(),
                                                                                   &temperature,
                                                                                   nullptr,
                                                                                   &frequency);

        while (true) {
            auto magneticFlux = OpenMagnetics::MagneticField::get_magnetic_flux(operationPoint.get_mutable_excitations_per_winding()[0].get_magnetizing_current().value(),
                                                                                neededTotalReluctance,
                                                                                numberTurnsPrimary,
                                                                                frequency);
            auto magneticFluxDensity = OpenMagnetics::MagneticField::get_magnetic_flux_density(magneticFlux, effectiveArea, frequency);
            auto magneticFieldStrength = OpenMagnetics::MagneticField::get_magnetic_field_strength(magneticFluxDensity, currentInitialPermeability, frequency);

            modifiedInitialPermeability = initial_permeability.get_initial_permeability(core.get_functional_description().get_material(),
                                                                                        &temperature,
                                                                                        &magneticFieldStrength.get_processed().value().get_mutable_offset(),
                                                                                        &frequency);

            if (fabs(currentInitialPermeability - modifiedInitialPermeability) < 1 || timeout == 0) {
                break;
            }
            else {
                currentInitialPermeability = modifiedInitialPermeability;
                timeout--;
            }
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
                    while(numberDistributedGaps > 3){
                        gappedCore = get_core_with_distributed_gapping(core, gapLength, numberDistributedGaps);
                        fringingFactorOneGap = reluctanceModel->get_gap_reluctance(gappedCore.get_gapping()[0])["fringing_factor"];
                        if (fringingFactorOneGap < constants.minimum_distributed_fringing_factor && numberDistributedGaps > 1) {
                            gapLength *= numberDistributedGaps;
                            numberDistributedGaps -= 2;
                            gapLength /= numberDistributedGaps;
                        }
                        else {
                            break;
                        }
                    }
                    while(true){
                        gappedCore = get_core_with_distributed_gapping(core, gapLength, numberDistributedGaps);
                        fringingFactorOneGap = reluctanceModel->get_gap_reluctance(gappedCore.get_gapping()[0])["fringing_factor"];
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

            if (fabs(neededTotalReluctance - reluctance) / neededTotalReluctance < 0.001 || timeout == 0 ) {
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
                if (increasingGap){
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

}
