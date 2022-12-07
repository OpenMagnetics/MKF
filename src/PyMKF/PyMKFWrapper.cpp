#include "pybind11_json.hpp"
#include <pybind11/pybind11.h>
#include "Constants.h"
#include "json.hpp"

using json = nlohmann::json;


#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;


py::dict get_constants(){
    py::dict dict;
    auto constants = OpenMagnetics::Constants();
    dict["residualGap"] = constants.residualGap;
    dict["minimumNonResidualGap"] = constants.minimumNonResidualGap;
    return dict;
}


PYBIND11_MODULE(PyMKF, m) {
    m.def("get_constants", &get_constants, "Returns the constants");
}
