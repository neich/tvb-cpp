//
// Created by imartin on 18-Oct-22.
//

#ifndef TVB_ROOT_CPP_PYTHON_API_H
#define TVB_ROOT_CPP_PYTHON_API_H

#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <pybind11/stl.h>

#include <tvb-cpp/simulator/simulate.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

void setWeights(py::EigenDRef<tvb::TArray2d> vref);
void setGlobalCoupling(tvb::Float g);
void setLengths(py::EigenDRef<tvb::TArray2d> vref, tvb::Float s);
void setIntegratorES(tvb::Float dt, py::EigenDRef<tvb::TArray1d> sigmas);
void setModel(std::string name);
void setModelParameter(std::string name, tvb::Float value);
void setModelParameter(std::string name, py::EigenDRef<tvb::TArray1d> value);
void setModelParameterSweep(std::string name, tvb::Float v_start, tvb::Float v_end, int n);
void printModelParameters();
void addRawMonitor(tvb::Float period, std::vector<int> voi);
void addTemporalAverageMonitor(tvb::Float period, std::vector<int> voi);
std::vector<std::tuple<py::array_t<tvb::Float>, py::array_t<tvb::Float>>> run_sim(tvb::Float t_start, tvb::Float t_end);

#endif //TVB_ROOT_CPP_PYTHON_API_H
