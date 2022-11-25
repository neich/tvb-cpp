#include <pybind11/pybind11.h>

#include "api.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

int add(int i, int j) {
    return i + j;
}

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: tvbcpp

        .. autosummary::
           :toctree: _generate

           add
           subtract
           model
    )pbdoc";

    m.def("set_weights", &setWeights, R"pbdoc(
    )pbdoc");

    m.def("set_lenghts", &setLengths, R"pbdoc(
    )pbdoc");

    m.def("set_integrator_es", &setIntegratorES, R"pbdoc(
    )pbdoc");

    m.def("set_model", &setModel, R"pbdoc(
    )pbdoc");

    m.def("set_model_parameter", static_cast<void (*)(std::string, tvb::Float)>(&setModelParameter), R"pbdoc(
    )pbdoc");

    m.def("set_model_parameter", static_cast<void (*)(std::string, py::EigenDRef<tvb::TArray1d>)>(&setModelParameter), R"pbdoc(
    )pbdoc");

    m.def("add_raw_monitor", &addRawMonitor, R"pbdoc(
    )pbdoc");

    m.def("add_average_monitor", &addAverageMonitor, R"pbdoc(
    )pbdoc");

    m.def("run_sim", &run_sim, R"pbdoc(
    )pbdoc");



#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
