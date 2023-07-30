//
// Created by imartin on 18-Oct-22.
//

#include "api.h"

#include <chrono>
#include <cassert>
#include <algorithm>

#include <tvb-root-cpp/simulator/integrators/euler_stochastic.h>
#include <tvb-root-cpp/simulator/noise.h>
#include <tvb-root-cpp/simulator/models/reduced_ww_ext.h>
#include <tvb-root-cpp/simulator/models/montbrio.h>
#include <tvb-root-cpp/simulator/models/zerlaut.h>
#include <tvb-root-cpp/simulator/noise.h>
#include <tvb-root-cpp/simulator/noise.h>
#include <tvb-root-cpp/simulator/simulator.h>


tvb::TArray2d weights;
tvb::TArray2d lengths;
float speed;
tvb::Integrator* integrator;
std::vector<tvb::Monitor*> monitors;
tvb::Model *model;
tvb::Monitor* monitor;
float dt = 0.1;
int N = 0;

void setWeights(py::EigenDRef<tvb::TArray2d> vref) {
    weights = vref;
    N = weights.rows();
    assert(("Matrix must be square!", N == weights.cols()));
}

void setLengths(py::EigenDRef<tvb::TArray2d> vref, float s) {
    lengths = vref;
    speed = s;
}

void setIntegratorES(float d, py::EigenDRef<tvb::TArray1d> sigmas) {
    integrator = new tvb::EulerStochastic(d, new tvb::Additive(sigmas, d));
    dt = d;
}

void setModel(std::string name) {
    assert(("Unknown number of regions, configure weight matrix first", N > 0));
    if (name == "ReducedWongWangExcInh")
        model = new tvb::ReducedWongWangExcInh(N);
    else if (name == "Montbrio")
        model = new tvb::Montbrio(weights.rows());
    else if (name == "ZerlautAdaptationFirstOrder")
        model = new tvb::ZerlautAdaptationFirstOrder(weights.rows());
    else if (name == "ZerlautAdptationSecondOrder")
        model = new tvb::ZerlautAdaptationSecondOrder(weights.rows());
    else throw std::runtime_error(string_format("Model wit name <%s> does not exist", name.c_str()));
}

void setModelParameter(std::string name, tvb::Float value) {
    model->set_param(name, value);
}

void setModelParameter(std::string name, py::EigenDRef<tvb::TArray1d> value) {
    model->set_param(name, value);
}

void printModelParameters() {
    std::cout << "Parameters:\n";
    for (auto pname: model->get_param_list())
        std::cout << pname << ": " << model->get_param_value(pname);
}

void addRawMonitor(float period, std::vector<int> voi) {
    monitors.push_back(new tvb::RawSubSample(period, dt, voi));
}

void addAverageMonitor(float period, std::vector<int> voi) {
    monitors.push_back(new tvb::TemporalAverage(weights.cols(), period, dt, voi));
}

py::array_t<tvb::Float> run_sim(float t_start, float t_end) {
    if (weights.rows() == 0)
        throw std::runtime_error("Weights matrix not initialized");

    if (weights.rows() != weights.cols())
        throw std::runtime_error("Weights matrix not square");

    if (lengths.rows() == 0)
        throw std::runtime_error("Lengths matrix not initialized");

    if (lengths.rows() != lengths.cols())
        throw std::runtime_error("Lengths matrix not square");

    if (model == nullptr)
        throw std::runtime_error("Model not initialized");

    if (integrator == nullptr)
        throw std::runtime_error("Integrator not initialized");

    tvb::Connectivity con(weights, lengths, speed);

    auto *coupling = new tvb::CouplingLinearSparse(con.weights(), con.delays(), model->cvars());

    tvb::Simulator simulator;

    int index = 0;
    std::vector<int> vois(model->state_vars().size());
    std::generate_n(vois.begin(), model->state_vars().size(), [&index]() { return index++;});
    monitor = new tvb::Raw(dt, vois);
    monitors.push_back(monitor);

    py::print("Starting simulation, t_start = ", t_start, ", t_end = ", t_end);
    auto start = std::chrono::high_resolution_clock::now();
    simulator.run(model, &con, integrator, monitors, coupling, t_start, t_end, nullptr);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start);
    py::print("Simulation ended in ", to_string(duration.count()), " msec");


    int n_records = monitor->getRecords().size();
    int n_voi = monitor->getRecords()[0].record.cols();
    int n_regions = monitor->getRecords()[0].record.rows();

    for (int rec = 0; rec < n_records; ++rec)
        if (!monitor->getRecords()[rec].record.allFinite()) {
            std::cerr << "Error in record " << rec << "\n";
            for (int voi = 0; voi < n_voi; ++voi)
                for (int reg = 0; reg < n_regions; ++reg)
                    std::cerr << monitor->getRecords()[rec].record(reg, voi) << ", ";
            std::cerr << std::endl;
            break;
        }

    size_t sizef = sizeof(tvb::Float);
    size_t size = n_records * n_voi * n_regions;
    size_t rec_size = n_voi * n_regions;
    auto *data = new tvb::Float[size];
    for (int rec = 0; rec < n_records; ++rec) {
        const tvb::TArray2d& record = monitor->getRecords()[rec].record;
        memcpy(&data[rec * rec_size], record.data(), rec_size*sizef);
    }

    py::capsule free_when_done(data, [](void *f) {
        auto *d = reinterpret_cast<tvb::Float *>(f);
        delete[] d;
    });

    return py::array_t<tvb::Float>(
            {n_records, n_voi, n_regions}, // shape
            {n_regions*n_voi*sizef, n_regions*sizef, sizef}, // C-style contiguous strides for double
            data, // the data pointer
            free_when_done);
}



