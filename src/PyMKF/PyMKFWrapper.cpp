#include "pybind11_json.hpp"
#include <pybind11/pybind11.h>
#include "Constants.h"
#include "Core.h"
#include "json.hpp"

using json = nlohmann::json;


#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

json get_core_data(json coreData){
    OpenMagnetics::Core core(coreData);
    return core;
}

py::dict get_constants(){
    py::dict dict;
    auto constants = OpenMagnetics::Constants();
    dict["residualGap"] = constants.residualGap;
    dict["minimumNonResidualGap"] = constants.minimumNonResidualGap;
    return dict;
}


PYBIND11_MODULE(PyMKF, m) {
    m.def("get_constants", &get_constants, "Returns the constants");
    m.def("get_core_data", &get_core_data, "Returns the processed data from a core");
}
