#include "pybind11_json.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "Constants.h"
#include "CoreWrapper.h"
#include "Reluctance.h"
#include "MagnetizingInductance.h"
#include "Utils.h"
#include "json.hpp"

using json = nlohmann::json;


#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

json get_core_data(json coreData, bool includeMaterialData=false){
    OpenMagnetics::CoreWrapper core(coreData, includeMaterialData);
    return core;
}

json get_material_data(json materialName){
    auto materialData = OpenMagnetics::find_data_by_name<OpenMagnetics::CoreMaterial>(materialName);
    return materialData;
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
                                                    json operationPointData){
    OpenMagnetics::CoreWrapper core(coreData);
    OpenMagnetics::WindingWrapper winding(windingData);
    OpenMagnetics::OperationPoint operationPoint(operationPointData);

    OpenMagnetics::MagnetizingInductance magnetizing_inductance(std::map<std::string, std::string>({{"reluctance", "ZHANG"}}));
    double magnetizingInductance = magnetizing_inductance.get_inductance_from_number_turns_and_gapping(core, winding, operationPointData);

    return magnetizingInductance;
}

// std::vector<CoreGap> get_gapping_from_number_turns_and_inductance(json core,
//                                                                   json winding,
//                                                                   json inputs,
//                                                                   std::string gappingTypeString, 
//                                                                   size_t decimals = 4);

// double get_number_turns_from_gapping_and_inductance(json core,
//                                                     json inputs);

PYBIND11_MODULE(PyMKF, m) {
    m.def("get_constants", &get_constants, "Returns the constants");
    m.def("get_material_data", &get_material_data, "Returns the all data about a given core material");
    m.def("get_core_data", &get_core_data, "Returns the processed data from a core");
    m.def("get_gap_reluctance", &get_gap_reluctance, "Returns the reluctance and fringing flux factor of a gap");
    m.def("get_gap_reluctance_model_information", &get_gap_reluctance_model_information, "Returns the information and average error for gap reluctance models");
    m.def("get_inductance_from_number_turns_and_gapping", &get_inductance_from_number_turns_and_gapping, "Returns the inductance of a core, given its number of turns and gapping");
}
