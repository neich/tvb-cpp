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

    m.def("set_global_coupling", &setGlobalCoupling, R"pbdoc(
    )pbdoc");

    m.def("set_lenghts", &setLengths, R"pbdoc(
    )pbdoc");

    m.def("set_integrator", &setIntegrator, R"pbdoc(
    )pbdoc");

    m.def("set_model", &setModel, R"pbdoc(
    )pbdoc");

    m.def("set_model_parameter", static_cast<void (*)(const std::string&, tvb::Float)>(&setModelParameter), R"pbdoc(
    )pbdoc");

    m.def("set_model_parameter", static_cast<void (*)(const std::string&, const py::EigenDRef<tvb::TArray1d>& value)>(&setModelParameter), R"pbdoc(
    )pbdoc");

    m.def("set_model_parameter_sweep", &setModelParameterSweep, R"pbdoc(
    )pbdoc");

    m.def("set_num_threads", &setNumThreads, R"pbdoc(
    )pbdoc");

    m.def("get_model_parameters", &getModelParameters, R"pbdoc(
    )pbdoc");

    m.def("add_raw_monitor", &addRawMonitor, R"pbdoc(
    )pbdoc");

    m.def("add_temporal_average_monitor", &addTemporalAverageMonitor, R"pbdoc(
    )pbdoc");

    m.def("add_bold_monitor", &addBOLDMonitor, R"pbdoc(
    )pbdoc");

    m.def("run_sim", &run_sim, R"pbdoc(
    )pbdoc");

    m.def("run_sweep", &run_sweep, R"pbdoc(
    )pbdoc");

    m.def("set_initial_state", &setInitialState, R"pbdoc(
    )pbdoc");



#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
