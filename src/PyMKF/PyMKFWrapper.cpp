#include "pybind11_json.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "Constants.h"
#include "CoreWrapper.h"
#include "Reluctance.h"
#include "Utils.h"
#include "json.hpp"

using json = nlohmann::json;


#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

json get_core_data(json coreData){
    OpenMagnetics::CoreWrapper core(coreData);
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


PYBIND11_MODULE(PyMKF, m) {
    m.def("get_constants", &get_constants, "Returns the constants");
    m.def("get_material_data", &get_material_data, "Returns the all data about a given core material");
    m.def("get_core_data", &get_core_data, "Returns the processed data from a core");
    m.def("get_gap_reluctance", &get_gap_reluctance, "Returns the reluctance and fringing flux factor of a gap");
    m.def("get_gap_reluctance_model_information", &get_gap_reluctance_model_information, "Returns the information and average error for gap reluctance models");
}
