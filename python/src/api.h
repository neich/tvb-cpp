//
// Created by imartin on 18-Oct-22.
//

#ifndef TVB_ROOT_CPP_PYTHON_API_H
#define TVB_ROOT_CPP_PYTHON_API_H

#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <pybind11/stl.h>

#include <tvb-root-cpp/simulator/simulate.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

void setWeights(py::EigenDRef<tvb::TArray2d> vref);
void setLengths(py::EigenDRef<tvb::TArray2d> vref, float s);
void setIntegratorES(float dt, py::EigenDRef<tvb::TArray1d> sigmas);
void setModel(std::string name);
void setModelParameter(std::string name, tvb::Float value);
void setModelParameter(std::string name, py::EigenDRef<tvb::TArray1d> value);
void addRawMonitor(float period, std::vector<int> voi);
void addAverageMonitor(float period, std::vector<int> voi);
py::array_t<tvb::Float> run_sim(float t_start, float t_end);

#endif //TVB_ROOT_CPP_PYTHON_API_H
