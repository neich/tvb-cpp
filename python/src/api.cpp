//
// Created by imartin on 18-Oct-22.
//

#include "api.h"

#include <chrono>
#include <cassert>
#include <algorithm>

#include <tvb-cpp/simulator/integrators/euler_stochastic.h>
#include <tvb-cpp/simulator/noise.h>
#include <tvb-cpp/simulator/models/reduced_ww_ext.h>
#include <tvb-cpp/simulator/models/montbrio.h>
#include <tvb-cpp/simulator/models/zerlaut.h>
#include <tvb-cpp/simulator/noise.h>
#include <tvb-cpp/simulator/noise.h>
#include <tvb-cpp/simulator/simulator.h>


tvb::TArray2d weights;
tvb::TArray2d lengths;
tvb::Float speed;
tvb::Integrator *integrator;
std::vector<tvb::Monitor *> monitors;
tvb::Model *model;
tvb::Monitor *monitor;
tvb::Float dt = 0.1;
int N = 0;
tvb::Float G = 1.0;

struct ParamSweep {
    tvb::Float v_start;
    tvb::Float v_end;
    int n;
};

std::unordered_map<std::string, ParamSweep> params;

void setWeights(py::EigenDRef<tvb::TArray2d> vref) {
    weights = vref;
    N = weights.rows();
    assert(("Matrix must be square!", N == weights.cols()));
}

void setGlobalCoupling(tvb::Float g) {
    assert(("G has to be a positive number!", g >= 0.0));
    G = g;
}

void setLengths(py::EigenDRef<tvb::TArray2d> vref, tvb::Float s) {
    lengths = vref;
    speed = s;
}

void setIntegratorES(tvb::Float d, py::EigenDRef<tvb::TArray1d> sigmas) {
    integrator = new tvb::EulerStochastic(d, new tvb::Additive(sigmas, d));
    dt = d;
}

void setModel(std::string name) {
    assert(("Unknown number of regions, configure weight matrix first", N > 0));
    if (name == "ReducedWongWangExcInh")
        model = new tvb::ReducedWongWangExcInh(N);
    else if (name == "Montbrio")
        model = new tvb::Montbrio(N);
    else if (name == "ZerlautAdaptationFirstOrder")
        model = new tvb::ZerlautAdaptationFirstOrder(N);
    else if (name == "ZerlautAdptationSecondOrder")
        model = new tvb::ZerlautAdaptationSecondOrder(N);
    else throw std::runtime_error(string_format("Model wit name <%s> does not exist", name.c_str()));

    model->configure();
}

void setModelParameter(std::string name, tvb::Float value) {
    model->set_param(name, value);
}

void setModelParameter(std::string name, py::EigenDRef<tvb::TArray1d> value) {
    model->set_param(name, value);
}

void setModelParameterSweep(std::string name, tvb::Float v_start, tvb::Float v_end, int n) {
    params[name] = ParamSweep(v_start, v_end, n);
}

void printModelParameters() {
    std::cout << "Parameters:\n";
    for (auto pname: model->get_param_list())
        std::cout << pname << ": " << model->get_param_value(pname);
}

void addRawMonitor(tvb::Float period, std::vector<int> voi) {
    monitors.push_back(new tvb::RawSubSample(period, dt, voi));
}

void addTemporalAverageMonitor(tvb::Float period, std::vector<int> voi) {
    monitors.push_back(new tvb::TemporalAverage(weights.cols(), period, dt, voi));
}

std::vector<std::tuple<py::array_t<tvb::Float>, py::array_t<tvb::Float>>>
run_sim(tvb::Float t_start, tvb::Float t_end) {
    if (weights.rows() == 0)
        throw std::runtime_error("Weights matrix not initialized");

    if (weights.rows() != weights.cols())
        throw std::runtime_error("Weights matrix not square");

    if (lengths.rows() == 0)
        lengths = tvb::TArray2d::Zero(weights.rows(), weights.cols());

    if (lengths.rows() != lengths.cols())
        throw std::runtime_error("Lengths matrix not square");

    if (model == nullptr)
        throw std::runtime_error("Model not initialized");

    if (integrator == nullptr)
        throw std::runtime_error("Integrator not initialized");

    model->init_dependant();

    tvb::Connectivity con(weights, lengths, speed);

    auto *coupling = new tvb::CouplingLinearSparse(con.weights(), con.delays(), model->cvars());
    coupling->setScale(G);

    tvb::Simulator simulator{};

    int index = 0;
    std::vector<int> vois(model->state_vars().size());
    std::generate_n(vois.begin(), model->state_vars().size(), [&index]() { return index++; });

    py::print("Starting simulation, t_start = ", t_start, ", t_end = ", t_end);
    auto start = std::chrono::high_resolution_clock::now();
    simulator.run(model, &con, integrator, monitors, coupling, t_start, t_end, nullptr);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start);
    py::print("Simulation ended in ", to_string(duration.count()), " msec");

    std::vector<std::tuple<py::array_t<tvb::Float>, py::array_t<tvb::Float>>> result;

    for (auto &monitor: monitors) {
        int n_records = monitor->getRecords().size();
        int n_voi = monitor->getRecords()[0].record.cols();
        int n_regions = monitor->getRecords()[0].record.rows();

//    for (int rec = 0; rec < n_records; ++rec)
//        if (!monitor->getRecords()[rec].record.allFinite()) {
//            std::cerr << "Error in record " << rec << "\n";
//            for (int voi = 0; voi < n_voi; ++voi)
//                for (int reg = 0; reg < n_regions; ++reg)
//                    std::cerr << monitor->getRecords()[rec].record(reg, voi) << ", ";
//            std::cerr << std::endl;
//            break;
//        }

        size_t sizef = sizeof(tvb::Float);
        size_t size = n_records * n_voi * n_regions;
        size_t rec_size = n_voi * n_regions;
        auto *data_monitor = new tvb::Float[size];
        auto *data_time = new tvb::Float[n_records];
        for (int rec = 0; rec < n_records; ++rec) {
            data_time[rec] = monitor->getRecords()[rec].time;
            const tvb::TArray2d &record = monitor->getRecords()[rec].record;
            memcpy(&data_monitor[rec * rec_size], record.data(), rec_size * sizef);
        }

        py::capsule free_when_done_monitor(data_monitor, [](void *f) {
            auto *d = reinterpret_cast<tvb::Float *>(f);
            delete[] d;
        });

        py::capsule free_when_done_time(data_time, [](void *f) {
            auto *d = reinterpret_cast<tvb::Float *>(f);
            delete[] d;
        });

        result.emplace_back(py::array_t<tvb::Float>(
                                    {n_records}, // shape
                                    {sizef}, // C-style contiguous strides for double
                                    data_time, // the data_monitor pointer
                                    free_when_done_time),
                            py::array_t<tvb::Float>(
                                    {n_records, n_voi, n_regions}, // shape
                                    {n_regions * n_voi * sizef, n_regions * sizef,
                                     sizef}, // C-style contiguous strides for double
                                    data_monitor, // the data_monitor pointer
                                    free_when_done_monitor));
    }

    return result;
}



