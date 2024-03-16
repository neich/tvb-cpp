//
// Created by imartin on 18-Oct-22.
//

#ifndef TVB_ROOT_CPP_PYTHON_API_H
#define TVB_ROOT_CPP_PYTHON_API_H

#include <pybind11/numpy.h>

#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <pybind11/stl.h>

#include <tvb-cpp/definitions.h>
#include <tvb-cpp/simulator/bold/BOLDModel.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

typedef py::array_t<tvb::Float> PyArray;
typedef std::tuple<PyArray, PyArray> PyTupleArray;
typedef typename std::vector<PyTupleArray> SimResult;
typedef typename std::tuple<std::string, tvb::Float> ParamScalar;
typedef typename std::tuple<std::string, tvb::TArray2d> ParamArray;
typedef typename std::vector<ParamScalar> ParamSet;
typedef typename std::vector<std::tuple<ParamSet, SimResult >> SweepResult;

const std::string field_prefix = "m_";


class BOLDModelWrapper;

void setWeights(py::EigenDRef<tvb::TArray2d> vref);
void setGlobalCoupling(tvb::Float g);
void setLengths(py::EigenDRef<tvb::TArray2d> vref, tvb::Float s);
void setIntegrator(const std::string &name, tvb::Float dt, py::EigenDRef<tvb::TArray1d> sigmas);
void setModel(const std::string &name);
void setModelParameter(const std::string& name, tvb::Float value);
void setModelParameter(const std::string& name, const py::EigenDRef<tvb::TArray1d>& value);
void setModelParameterSweep(const std::string& name, tvb::Float v_start, tvb::Float v_end, int n);
void setInitialState(py::EigenDRef<tvb::TArray2d> is);
void setNumThreads(int n);
std::vector<std::tuple<std::string, tvb::TArray1d>> getModelParameters();
void addRawMonitor(tvb::Float period, std::vector<int> voi);
void addTemporalAverageMonitor(tvb::Float period, std::vector<int> voi);
BOLDModelWrapper create_bold(const std::string& type, tvb::Float tr);
SimResult run_sim(tvb::Float t_start, tvb::Float t_end);
SweepResult run_sweep(tvb::Float t_start, tvb::Float t_end);


class BOLDModelWrapper {
    tvb::BOLDModel *m_bm;

public:
    BOLDModelWrapper(tvb::BOLDModel *bm): m_bm(bm) {}

    void set_param(const std::string& pname, tvb::Float value)  {
        m_bm->set_param(field_prefix + pname, value);
    }

    [[nodiscard]] std::vector<std::string> get_param_list() const {
        std::vector<std::string> plist = m_bm->get_param_list();
        std::vector<std::string> result;
        std::string prefix = "m_";
        std::transform(plist.begin(), plist.end(), std::back_inserter(result), [&prefix](const std::string& p) { return p.substr(field_prefix.size()); });
        return result;
    }

    tvb::Float get_param(const std::string& param) const {
        return m_bm->get_param(field_prefix + param);
    }

    [[nodiscard]] PyTupleArray compute_bold(const tvb::TArray2d& ts, tvb::Float ts_dt) const;
};
#endif //TVB_ROOT_CPP_PYTHON_API_H
