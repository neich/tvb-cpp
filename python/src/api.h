//
// Created by imartin on 18-Oct-22.
//

#ifndef TVB_ROOT_CPP_PYTHON_API_H
#define TVB_ROOT_CPP_PYTHON_API_H

#include <pybind11/numpy.h>

#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <pybind11/stl.h>

#include <tvb-cpp/definitions.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

typedef typename std::vector<std::tuple<py::array_t<tvb::Float>, py::array_t<tvb::Float>>> SimResult;
typedef typename std::tuple<std::string, tvb::Float> ParamScalar;
typedef typename std::tuple<std::string, tvb::TArray2d> ParamArray;
typedef typename std::vector<ParamScalar> ParamSet;
typedef typename std::vector<std::tuple<ParamSet, SimResult >> SweepResult;

void setWeights(py::EigenDRef<tvb::TArray2d> vref);
void setGlobalCoupling(tvb::Float g);
void setLengths(py::EigenDRef<tvb::TArray2d> vref, tvb::Float s);
void setIntegratorES(tvb::Float dt, py::EigenDRef<tvb::TArray1d> sigmas);
void setModel(std::string name);
void setModelParameter(const std::string& name, tvb::Float value);
void setModelParameter(const std::string& name, const py::EigenDRef<tvb::TArray1d>& value);
void setModelParameterSweep(const std::string& name, tvb::Float v_start, tvb::Float v_end, int n);
void setNumThreads(int n);
void printModelParameters();
void addRawMonitor(tvb::Float period, std::vector<int> voi);
void addTemporalAverageMonitor(tvb::Float period, std::vector<int> voi);
SimResult run_sim(tvb::Float t_start, tvb::Float t_end);
SweepResult run_sweep(tvb::Float t_start, tvb::Float t_end);

#endif //TVB_ROOT_CPP_PYTHON_API_H
