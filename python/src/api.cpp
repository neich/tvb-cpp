//
// Created by imartin on 18-Oct-22.
//

#include "api.h"

#include <chrono>
#include <cassert>
#include <algorithm>

#include <tvb-cpp/simulator/integrators/euler_stochastic.h>
#include <tvb-cpp/simulator/integrators/euler_deterministic.h>
#include <tvb-cpp/simulator/models/reduced_ww_ext.h>
#include <tvb-cpp/simulator/models/montbrio.h>
#include <tvb-cpp/simulator/models/zerlaut.h>
#include <tvb-cpp/simulator/noise.h>
#include <tvb-cpp/simulator/noise.h>
#include <tvb-cpp/simulator/simulator.h>
#include <tvb-cpp/tools/threadpool.h>
#include <tvb-cpp/simulator/monitors/bold_tvb.h>

#ifdef __unix__
# include <unistd.h>

#elif defined _WIN32
# include <windows.h>
#define sleep(x) Sleep(1000 * (x))
#endif

using namespace std;

tvb::TArray1d EMPTY(0);

tvb::TArray2d weights;
tvb::TArray2d lengths;
tvb::Float speed;
string integrator_name;
tvb::TArray1d int_sigmas = EMPTY;
tvb::Float int_dt = 0.1;
std::vector<tvb::Monitor *> monitors;
std::string model_name{};
tvb::Monitor *monitor;
tvb::Float dt = 0.1;
int N = 0;
tvb::Float G = 1.0;
int num_threads = 1;

struct ParamSweep {
    tvb::Float v_start;
    tvb::Float v_end;
    int n;
};


void checkState();

SimResult genSimResultFromMonitors(const std::vector<tvb::Monitor *> &mntrs);

tvb::Model *genModel(const string &name);

std::unordered_map<std::string, ParamSweep> params_sweep;
std::vector<ParamScalar> params_scalar;
std::vector<ParamArray> params_array;

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

void setIntegrator(const string& name, tvb::Float d, py::EigenDRef<tvb::TArray1d> sgms = EMPTY) {
    integrator_name = name;
    int_sigmas = sgms;
    int_dt = d;
}

void setModel(const std::string& name) {
    model_name = name;
}

tvb::Model *genModel(const string &name) {
    assert(("Unknown number of regions, configure weight matrix first", N > 0));
    tvb::Model *model;
    if (name == "ReducedWongWangExcInh")
        model = new tvb::ReducedWongWangExcInh(N);
    else if (name == "Montbrio")
        model = new tvb::Montbrio(N);
    else if (name == "ZerlautAdaptationFirstOrder")
        model = new tvb::ZerlautAdaptationFirstOrder(N);
    else if (name == "ZerlautAdptationSecondOrder")
        model = new tvb::ZerlautAdaptationSecondOrder(N);
    else throw runtime_error(string_format("Model with name <%s> does not exist", name.c_str()));

    model->configure();

    return model;
}
tvb::Integrator *genIntegrator(const string &name) {
    tvb::Integrator *integrator_ret;
    if (name == "EulerDeterministic")
        integrator_ret = new tvb::EulerDeterministic(N);
    else if (name == "EulerStochastic")
        integrator_ret = new tvb::EulerStochastic(int_dt, new tvb::Additive(int_sigmas, int_dt));
    else throw runtime_error(string_format("Integrator with name <%s> does not exist", name.c_str()));

    return integrator_ret;
}

void setModelParameter(const std::string& name, tvb::Float value) {
    params_scalar.emplace_back(name, value);
}

void setModelParameter(const std::string& name, const py::EigenDRef<tvb::TArray1d>& value) {
    params_array.emplace_back(name, value);
}

void setModelParameterSweep(const std::string& name, tvb::Float v_start, tvb::Float v_end, int n) {
    params_sweep[name] = ParamSweep(v_start, v_end, n);
}

void setNumThreads(int n) {
    assert(("It has to be a number greater than 0!", n > 0));
    num_threads = n;
}

void printModelParameters() {
    std::cout << "Parameters:\n";
}

void addRawMonitor(tvb::Float period, std::vector<int> voi) {
    monitors.push_back(new tvb::RawSubSample(period, dt, voi));
}

void addTemporalAverageMonitor(tvb::Float period, std::vector<int> voi) {
    monitors.push_back(new tvb::TemporalAverage(weights.cols(), period, dt, voi));
}

void addBOLDMonitor(tvb::Float period, std::vector<int> voi) {
    monitors.push_back(new tvb::BoldTVB(weights.cols(), period, dt, voi));
}

void simulate(const tvb::Model *model,
              const tvb::Connectivity *con,
              const tvb::Integrator *integrator,
              const std::vector<tvb::Monitor*> &monitors,
              tvb::Coupling *coupling,
              tvb::Float t_start,
              tvb::Float t_end) {

    tvb::Simulator simulator{};
    simulator.run(model, con, integrator, monitors, coupling, t_start, t_end, nullptr);
}


SimResult
run_sim(tvb::Float t_start, tvb::Float t_end) {
    if (!params_sweep.empty())
        throw std::runtime_error("Cannot run a single simulation, a parameter sweep has been defined!");

    checkState();

    tvb::Model *model = genModel(model_name);
    for (auto const &p: params_scalar)
        model->set_param(std::get<0>(p), std::get<1>(p));
    for (auto const &p: params_array)
        model->set_param(std::get<0>(p), std::get<1>(p));
    model->init_dependant();

    tvb::Connectivity con(weights, lengths, speed);

    auto *coupling = new tvb::CouplingLinearSparse(con.weights(), con.delays(), model->cvars());
    coupling->setScale(G);

    tvb::Simulator simulator{};

//    int index = 0;
//    std::vector<int> vois(model->state_vars().size());
//    std::generate_n(vois.begin(), model->state_vars().size(), [&index]() { return index++; });

    py::print("Starting simulation, t_start = ", t_start, ", t_end = ", t_end);
    auto start = std::chrono::high_resolution_clock::now();
    auto *integrator = genIntegrator(integrator_name);
    simulate(model, &con, integrator, monitors, coupling, t_start, t_end);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start);
    py::print("Simulation ended in ", to_string(duration.count()), " msec");

    SimResult result = genSimResultFromMonitors(monitors);

    return result;
}

SimResult genSimResultFromMonitors(const std::vector<tvb::Monitor *> &mntrs) {
    SimResult result;

    for (auto &m: mntrs) {
        int n_records = m->getRecords().size();
        if (n_records == 0) {
            py::print("Monitor empty!");
            continue;
        }

        int n_voi = m->getRecords()[0].record.cols();
        int n_regions = m->getRecords()[0].record.rows();

        size_t sizef = sizeof(tvb::Float);
        size_t size = n_records * n_voi * n_regions;
        size_t rec_size = n_voi * n_regions;
        auto *data_monitor = new tvb::Float[size];
        auto *data_time = new tvb::Float[n_records];
        for (int rec = 0; rec < n_records; ++rec) {
            data_time[rec] = m->getRecords()[rec].time;
            const tvb::TArray2d &record = m->getRecords()[rec].record;
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

        result.emplace_back(pybind11::array_t<tvb::Float>(
                                    {n_records}, // shape
                                    {sizef}, // C-style contiguous strides for double
                                    data_time, // the data_monitor pointer
                                    free_when_done_time),
                            pybind11::array_t<tvb::Float>(
                                    {n_records, n_voi, n_regions}, // shape
                                    {n_regions * n_voi * sizef, n_regions * sizef,
                                     sizef}, // C-style contiguous strides for double
                                    data_monitor, // the data_monitor pointer
                                    free_when_done_monitor));
    }
    return result;
}

void checkState() {
    if (weights.rows() == 0)
        throw runtime_error("Weights matrix not initialized");

    if (weights.rows() != weights.cols())
        throw runtime_error("Weights matrix not square");

    if (lengths.rows() == 0)
        lengths = tvb::TArray2d::Zero(weights.rows(), weights.cols());

    if (lengths.rows() != lengths.cols())
        throw runtime_error("Lengths matrix not square");

    if (model_name.empty())
        throw runtime_error("Model not initialized");

    if (integrator_name.empty())
        throw runtime_error("Integrator not initialized");
}

SweepResult run_sweep(tvb::Float t_start, tvb::Float t_end) {
    assert(("No sweep defined!", !params_sweep.empty()));

    try {
        checkState();

        std::vector<ParamSet> param_combs(1);
        for (auto const &p: params_sweep) {
            std::vector<ParamSet> new_param_combs;
            for (auto v: tvb::range(p.second.v_start, p.second.v_end, p.second.n)) {
                for (auto const &pc: param_combs) {
                    new_param_combs.push_back(pc);
                    new_param_combs.back().emplace_back(p.first, v);
                }
            }
            param_combs = new_param_combs;
        }
        for (auto &pc: param_combs) {
            for (auto const &p: params_scalar)
                pc.emplace_back(std::get<0>(p), std::get<1>(p));
        }


        tvb::Connectivity con(weights, lengths, speed);

        std::vector<std::vector<tvb::Monitor *>> sim_results;
        for (auto &pc: param_combs) {
            sim_results.emplace_back();
            for (auto const *m: monitors) {
                sim_results.back().push_back(m->clone());
            }
        }

        auto *integrator = genIntegrator(integrator_name);

        py::print(string_format("Creating thread pool with %d threads", num_threads));
        tvb::ThreadPool<int> tp(num_threads);
        tp.start();

        int nsim = 0;
        for (auto &pc: param_combs) {

            tvb::Model *model = genModel(model_name);
            for (auto const &p: pc)
                model->set_param(std::get<0>(p), std::get<1>(p));
            model->init_dependant();

            auto *coupling = new tvb::CouplingLinearSparse(con.weights(), con.delays(), model->cvars());
            coupling->setScale(G);

            std::string msg = string_format("Starting sweep (%d of %d) for:", nsim, param_combs.size());
            for (auto const &p: pc)
                msg += string_format(" %s=%f", std::get<0>(p).c_str(), std::get<1>(p));

            py::print(msg);

            const tvb::Connectivity *con_ref = &con;
            const std::vector<tvb::Monitor *> &sim_monitors = sim_results[nsim];
            tp.queue_job([model, con_ref, integrator, sim_monitors, coupling, t_start, t_end, nsim] {
                simulate(model, con_ref, integrator, sim_monitors, coupling, t_start, t_end);
                return nsim;
            });
            nsim++;
        }

        SweepResult return_values;
        while (!tp.finished()) {
            std::optional<int> op = tp.get_result();
            if (op.has_value()) {
                int nsim = op.value();
                return_values.emplace_back(param_combs[nsim], genSimResultFromMonitors(sim_results[nsim]));
            }
            sleep(1);
        }

        tp.stop();

        return return_values;
    } catch (std::runtime_error e) {
        py::print(string_format("TVB error:%s", e.what()));
    }
}