#include "pybind11_json.hpp"
#include <pybind11/pybind11.h>
#include "Constants.h"
#include "json.hpp"

using json = nlohmann::json;


#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

int add(int i, int j) {
    return i + j;
}

namespace py = pybind11;

PYBIND11_MODULE(PyMKF, m) {
    py::class_<OpenMagnetics::Constants>(m, "Constants")
        .def(py::init())
        .def_readonly("residual_gap", &OpenMagnetics::Constants::residualGap);
}
