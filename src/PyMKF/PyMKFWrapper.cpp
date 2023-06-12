#include "pybind11_json.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "pybind11_json/pybind11_json.hpp"
#include "Constants.h"
#include "Defaults.h"
#include <MAS.hpp>
#include "InputsWrapper.h"
#include "CoreWrapper.h"
#include "Reluctance.h"
#include "MagnetizingInductance.h"
#include "CoreLosses.h"
#include "CoreTemperature.h"
#include "Utils.h"
#include "json.hpp"
#include <pybind11/embed.h>
#include <magic_enum.hpp>

using json = nlohmann::json;


#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

json get_core_data(json coreData, bool includeMaterialData=false){
    OpenMagnetics::CoreWrapper core(coreData, includeMaterialData);
    return core;
}

json get_material_data(json materialName){
    auto materialData = OpenMagnetics::find_core_material_by_name(materialName);
    return materialData;
}

json get_shape_data(json shapeName){
    auto shapeData = OpenMagnetics::find_core_shape_by_name(shapeName);
    return shapeData;
}

json get_available_shape_families(){
    return magic_enum::enum_names<OpenMagnetics::CoreShapeFamily>();
}

json get_available_core_materials(){
    return OpenMagnetics::get_material_names();
}
json get_available_core_shapes(){
    return OpenMagnetics::get_shape_names();
}

py::dict get_constants(){
    py::dict dict;
    auto constants = OpenMagnetics::Constants();
    dict["residualGap"] = constants.residualGap;
    dict["minimumNonResidualGap"] = constants.minimumNonResidualGap;
    return dict;
}

json get_gap_reluctance(json coreGapData, std::string modelNameString){

    auto modelName = magic_enum::enum_cast<OpenMagnetics::ReluctanceModels>(modelNameString);

    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(modelName.value());
    OpenMagnetics::CoreGap coreGap(coreGapData);

    auto coreGapResult = reluctanceModel->get_gap_reluctance(coreGap);
    return coreGapResult;
}

json get_gap_reluctance_model_information(){
    py::dict dict;
    dict["information"] = OpenMagnetics::ReluctanceModel::get_models_information();
    dict["errors"] = OpenMagnetics::ReluctanceModel::get_models_errors();
    dict["internal_links"] = OpenMagnetics::ReluctanceModel::get_models_internal_links();
    dict["external_links"] = OpenMagnetics::ReluctanceModel::get_models_external_links();
    return dict;
}

double get_inductance_from_number_turns_and_gapping(json coreData,
                                                    json windingData,
                                                    json operationPointData,
                                                    json modelsData){
    OpenMagnetics::CoreWrapper core(coreData);
    OpenMagnetics::WindingWrapper winding(windingData);
    OpenMagnetics::OperationPoint operationPoint(operationPointData);

    std::map<std::string, std::string> models = modelsData.get<std::map<std::string, std::string>>();

    OpenMagnetics::MagnetizingInductance magnetizing_inductance(models);
    double magnetizingInductance = magnetizing_inductance.get_inductance_from_number_turns_and_gapping(core, winding, &operationPoint);

    return magnetizingInductance;
}


double get_number_turns_from_gapping_and_inductance(json coreData,
                                                    json inputsData,    
                                                    json modelsData){
    OpenMagnetics::CoreWrapper core(coreData);
    OpenMagnetics::InputsWrapper inputs(inputsData);

    std::map<std::string, std::string> models = modelsData.get<std::map<std::string, std::string>>();

    OpenMagnetics::MagnetizingInductance magnetizing_inductance(models);
    double numberTurns = magnetizing_inductance.get_number_turns_from_gapping_and_inductance(core, &inputs);

    return numberTurns;
}

py::list get_gapping_from_number_turns_and_inductance(json coreData,
                                                      json windingData,
                                                      json inputsData,
                                                      std::string gappingTypeString,
                                                      size_t decimals,
                                                      json modelsData){
    OpenMagnetics::CoreWrapper core(coreData);
    OpenMagnetics::WindingWrapper winding(windingData);
    OpenMagnetics::InputsWrapper inputs(inputsData);

    std::map<std::string, std::string> models = modelsData.get<std::map<std::string, std::string>>();
    OpenMagnetics::GappingType gappingType = magic_enum::enum_cast<OpenMagnetics::GappingType>(gappingTypeString).value();

    OpenMagnetics::MagnetizingInductance magnetizing_inductance(models);
    std::vector<OpenMagnetics::CoreGap> gapping = magnetizing_inductance.get_gapping_from_number_turns_and_inductance(core,
                                                                                                                      winding,
                                                                                                                      &inputs,
                                                                                                                      gappingType, 
                                                                                                                      decimals);

    py::list gappingAux;
    for (auto &gap: gapping) {
        json aux;
        to_json(aux, gap);
        py::dict pyAux = aux;
        gappingAux.append(pyAux);
    }
    return gappingAux;

}

json get_steinmetz_coefficients(std::string material, double frequency){
    auto steinmetz = py::module::import("steinmetz");
    json coefficients = steinmetz.attr("fit")(material, frequency).cast<std::map<std::string, double>>();


    return coefficients;
}

json get_core_losses(json coreData,
                     json windingData,
                     json inputsData,    
                     json modelsData){
    OpenMagnetics::CoreWrapper core(coreData);
    OpenMagnetics::WindingWrapper winding(windingData);
    OpenMagnetics::InputsWrapper inputs(inputsData);
    auto operationPoint = inputs.get_operation_point(0);
    OpenMagnetics::OperationPointExcitation excitation = operationPoint.get_excitations_per_winding()[0];
    double magnetizingInductance = OpenMagnetics::InputsWrapper::get_requirement_value(inputs.get_design_requirements().get_magnetizing_inductance());
    if (!excitation.get_current()) {
        auto magnetizingCurrent = OpenMagnetics::InputsWrapper::get_magnetizing_current(excitation, magnetizingInductance);
        excitation.set_current(magnetizingCurrent);
        operationPoint.get_mutable_excitations_per_winding()[0] = excitation;
    }

    auto defaults = OpenMagnetics::Defaults();

    bool enableTemperatureConvergence = false;

    std::map<std::string, std::string> models = modelsData.get<std::map<std::string, std::string>>();

    auto coreLossesModelName = defaults.coreLossesModelDefault;
    if (models.find("coreLosses") != models.end()) {
        coreLossesModelName = magic_enum::enum_cast<OpenMagnetics::CoreLossesModels>(models["coreLosses"]).value();
    }
    auto coreTemperatureModelName = defaults.coreTemperatureModelDefault;
    if (models.find("coreTemperature") != models.end()) {
        coreTemperatureModelName = magic_enum::enum_cast<OpenMagnetics::CoreTemperatureModels>(models["coreTemperature"]).value();
    }

    auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(models);
    auto coreTemperatureModel = OpenMagnetics::CoreTemperatureModel::factory(coreTemperatureModelName);

    OpenMagnetics::MagnetizingInductance magnetizing_inductance(models);

    double temperature = operationPoint.get_conditions().get_ambient_temperature();
    double temperatureAfterLosses = temperature;
    OpenMagnetics::ElectromagneticParameter magneticFluxDensity;
    json result;

    do {
        temperature = temperatureAfterLosses;

        excitation = operationPoint.get_excitations_per_winding()[0];
        operationPoint.get_mutable_conditions().set_ambient_temperature(temperature);

        magneticFluxDensity = magnetizing_inductance.get_inductance_and_magnetic_flux_density(core, winding, &operationPoint).second;
        excitation.set_magnetic_flux_density(magneticFluxDensity);

        result = coreLossesModel->get_core_losses(core, excitation, temperature);

        auto temperatureResult = coreTemperatureModel->get_core_temperature(core, result["totalLosses"].get<double>(), temperature);
        temperatureAfterLosses = temperatureResult["maximumTemperature"];
    } while (fabs(temperature - temperatureAfterLosses) / temperatureAfterLosses >= 0.01 && enableTemperatureConvergence);

    result["magneticFluxDensityPeak"] = magneticFluxDensity.get_processed().value().get_peak().value();
    result["magneticFluxDensityAcPeak"] = magneticFluxDensity.get_processed().value().get_peak().value() - magneticFluxDensity.get_processed().value().get_offset();
    result["voltageRms"] = operationPoint.get_mutable_excitations_per_winding()[0].get_voltage().value().get_processed().value().get_rms().value();
    result["currentRms"] = operationPoint.get_mutable_excitations_per_winding()[0].get_current().value().get_processed().value().get_rms().value();
    result["apparentPower"] = operationPoint.get_mutable_excitations_per_winding()[0].get_voltage().value().get_processed().value().get_rms().value() * operationPoint.get_mutable_excitations_per_winding()[0].get_current().value().get_processed().value().get_rms().value();
    result["maximumCoreTemperature"] = temperatureAfterLosses;
    result["maximumCoreTemperatureRise"] = temperatureAfterLosses - operationPoint.get_conditions().get_ambient_temperature();

    if (models["coreLosses"] == "ROSHEN") {
        result["_hysteresisMajorLoopTop"] = coreLossesModel->_hysteresisMajorLoopTop;
        result["_hysteresisMajorLoopBottom"] = coreLossesModel->_hysteresisMajorLoopBottom;
        result["_hysteresisMajorH"] = coreLossesModel->_hysteresisMajorH;
        result["_hysteresisMinorLoopTop"] = coreLossesModel->_hysteresisMinorLoopTop;
        result["_hysteresisMinorLoopBottom"] = coreLossesModel->_hysteresisMinorLoopBottom;
    }

    return result;
}


json get_core_losses_model_information(std::string material){
    py::dict dict;
    dict["information"] = OpenMagnetics::CoreLossesModel::get_models_information();
    dict["errors"] = OpenMagnetics::CoreLossesModel::get_models_errors();
    dict["internal_links"] = OpenMagnetics::CoreLossesModel::get_models_internal_links();
    dict["external_links"] = OpenMagnetics::CoreLossesModel::get_models_external_links();
    dict["available_models"] = OpenMagnetics::CoreLossesModel::get_methods(material);
    return dict;
}


json get_core_temperature_model_information(){
    py::dict dict;
    dict["information"] = OpenMagnetics::CoreTemperatureModel::get_models_information();
    dict["errors"] = OpenMagnetics::CoreTemperatureModel::get_models_errors();
    dict["internal_links"] = OpenMagnetics::CoreTemperatureModel::get_models_internal_links();
    dict["external_links"] = OpenMagnetics::CoreTemperatureModel::get_models_external_links();
    return dict;
}


PYBIND11_MODULE(PyMKF, m) {
    m.def("get_constants", &get_constants, "Returns the constants");
    m.def("get_available_shape_families", &get_available_shape_families, "Returns the available shape families");
    m.def("get_material_data", &get_material_data, "Returns the all data about a given core material");
    m.def("get_shape_data", &get_shape_data, "Returns the all data about a given core shape");
    m.def("get_available_core_materials", &get_available_core_materials, "Returns the names of all the core materials");
    m.def("get_available_core_shapes", &get_available_core_shapes, "Returns the names of all the core shapes");
    m.def("get_core_data", &get_core_data, "Returns the processed data from a core");
    m.def("get_gap_reluctance", &get_gap_reluctance, "Returns the reluctance and fringing flux factor of a gap");
    m.def("get_gap_reluctance_model_information", &get_gap_reluctance_model_information, "Returns the information and average error for gap reluctance models");
    m.def("get_inductance_from_number_turns_and_gapping", &get_inductance_from_number_turns_and_gapping, "Returns the inductance of a core, given its number of turns and gapping");
    m.def("get_number_turns_from_gapping_and_inductance", &get_number_turns_from_gapping_and_inductance, "Returns the number of turns needed to achieve a given inductance with a given gapping");
    m.def("get_gapping_from_number_turns_and_inductance", &get_gapping_from_number_turns_and_inductance, "Returns the gapping needed to achieve a given inductance with a given number of turns");
    m.def("get_steinmetz_coefficients", &get_steinmetz_coefficients, "");
    m.def("get_core_losses", &get_core_losses, "Returns the core losses according to given model");
    m.def("get_core_losses_model_information", &get_core_losses_model_information, "Returns the information and average error for core losses models");
    m.def("get_core_temperature_model_information", &get_core_temperature_model_information, "Returns the information and average error for core temperature models");
}
