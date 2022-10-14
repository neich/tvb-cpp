#include <pybind11/pybind11.h>

#include <tvb-root-cpp/simulator/models/reduced_ww_ext.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

int add(int i, int j) {
    return i + j;
}

void model(int N, float value) {
    auto *model = new tvb::ReducedWongWangExcInh(N);
    model->set_param_fill("G", value);
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

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");


    m.def("model", &model, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");


#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
