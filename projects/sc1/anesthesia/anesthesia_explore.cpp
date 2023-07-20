//    Copyright 2020-2021 Ignacio Martín <ignacio.martin@udg.edu>
//
//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at
//
//            http://www.apache.org/licenses/LICENSE-2.0
//
//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
//    limitations under the License.

#include <string>
#include <chrono>

#include <tvb-cpp/tools/npz_tools.h>
#include <tvb-cpp/simulator/simulate.h>
#include <tvb-cpp/simulator/monitor.h>
#include <tvb-cpp/simulator/models/reduced_ww_ext.h>
#include <tvb-cpp/simulator/models/montbrio.h>
#include <tvb-cpp/simulator/integrators/euler_stochastic.h>
#include <tvb-cpp/tools/threadpool.h>
#include "tvb-cpp/tools/csv_tools.h"
#include <tvb-cpp/tools/npy.h>
#include "tvb-cpp/simulator/integrators/euler_deterministic.h"
#include "zerlaut_gaba.h"
#include "tvb-cpp/simulator/simulator.h"
#include "tvb-cpp/simulator/monitors/bold_tvb.h"
#include "tvb-cpp/simulator/monitors/bold_BalloonWindkessel.h"
#include <tvb-cpp/tools/algo/fic/functions/balance_fic.h>

#include <chrono>
#include <filesystem>
#include <thread>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <utility>

using namespace boost::program_options;
using namespace std::filesystem;
using namespace tvb;
using namespace std::chrono;

float base_value = 0.3772258064;

struct SweepParam {
    string name;
    std::vector<float> values;
};

struct Parameter {
    string name;
    float value;

    Parameter(string name, float value) : name(std::move(name)), value(value) {}

    bool operator==(const std::string &v) const { return name == v; }
};

struct RunParams {
    float dt = 0.1;
    float t_start = 0.0;
    float t_end = 10000.0;
    bool force_output = false;
    int voi;
    float value_base;
    std::vector<float> sigmas{0.0, 0.0, 0.0, 0.0};
    tvb::Monitor *monitor = nullptr;
    std::vector<Parameter> params;
    string file_out;
    string file_prefix;
    string path_out;
    string file_weights;
    string file_lengths;
    string gaba_vector;
    float speed = 1e6;
    float G = 1.0;

    RunParams() = default;

    RunParams(std::vector<Parameter> params, tvb::Monitor *monitor) : params(std::move(params)), monitor(monitor) {}

    void init(const variables_map &vm) {
        this->voi = vm["var-of-interest"].as<int>();
        this->value_base = vm["value-base"].as<float>();
        this->file_weights = vm["sc-matrix"].as<std::string>();
        this->gaba_vector = vm["gaba-vector"].as<std::string>();
        this->t_start = vm["time-start"].as<float>();
        this->t_end = vm["time-end"].as<float>();
        this->dt = vm["dt"].as<float>();
        if (vm.count("length-matrix"))
            this->file_lengths = vm["length-matrix"].as<std::string>();
        // this->file_prefix = vm["out-file-prefix"].as<std::string>();
        this->path_out = vm["out-path"].as<std::string>();
        this->speed = vm["speed"].as<float>();
        this->force_output = vm["force-output"].as<bool>();
    }
};


//void save_fig(tvb::Monitor *monitor, const string &file_prefix) {
//    if (monitor == nullptr) return;
//
//    int n_records = monitor->getRecords().size();
//    int N = monitor->getRecords()[0].record.rows();
//    std::vector<std::vector<Float>> y_plot(N, std::vector<Float>(n_records));
//    for (unsigned t = 0; t < n_records; ++t) {
//        const Monitor::Record &r = monitor->getRecords()[t];
//        for (unsigned n = 0; n < N; ++n)
//            y_plot[n][t] = r.record(n, 0);
//    }
//
//    // tvb::csv_save("./test_simulationRWW_TVB_CPP.csv", y_plot);
//
//    // Plot line from given x and y data. Color is selected automatically.
//    std::vector<Float> ls(n_records);
//    std::transform(monitor->getRecords().begin(), monitor->getRecords().end(), ls.begin(),
//                   [](const Monitor::Record &r) { return r.time / 1000; });
//
//    {
//        std::unique_lock<std::mutex> lock(matplotlib_mutex);
//
//        for (unsigned n = 0; n < N; ++n) {
//            plt::plot(ls, y_plot[n]);
//        }
//        // Plot a red dashed line from given x and y data.
//        // plt::plot(x, w,"r--");
//        // Plot a line whose name will show up as "log(x)" in the legend.
//
//        string title = "Zerlaut GABA";
//
//        plt::title(title);
//        plt::ylabel("E");
//        plt::xlabel("Seconds");
////        plt::ylim(-0.01, 0.05);
////        plt::axhline(base_value - 0.005);
////        plt::axhline(base_value + 0.005);
//        // Save the image (file format is determined by the extension)
//
//        plt::save(file_prefix + ".png", 300);
//
//        plt::clf();
//    }
//}

void save_cvs(tvb::Monitor *monitor, const string &file_prefix) {
    if (monitor == nullptr) return;

    int n_records = monitor->getRecords().size();
    int N = monitor->getRecords()[0].record.rows();
    std::vector<std::vector<Float>> y_plot(N, std::vector<Float>(n_records));
    for (unsigned t = 0; t < n_records; ++t) {
        const Monitor::Record &r = monitor->getRecords()[t];
        for (unsigned n = 0; n < N; ++n)
            y_plot[n][t] = r.record(n, 0);
    }

    // Plot line from given x and y data. Color is selected automatically.
    std::vector<Float> ls(n_records);
    std::transform(monitor->getRecords().begin(), monitor->getRecords().end(), ls.begin(),
                   [](const Monitor::Record &r) { return r.time / 1000; });

    y_plot.insert(y_plot.begin(), ls);

    tvb::csv_save(file_prefix + ".cvs", y_plot, true);

}

void
load_data(const string &file_weights, const string &file_lengths, const string &file_gaba, TArray2d &C, TArray2d &tl,
          TArray1d &gaba_vector) {
    if (file_weights.ends_with(".csv"))
        C = tvb::csv_load(file_weights);
    else if (file_weights.ends_with(".npz"))
        C = tvb::npz2Matrixd(file_weights, "SC");
    else if (file_weights.ends_with(".npy"))
        C = tvb::npy2Matrixd(file_weights);
    else
        throw std::runtime_error(string_format("Unknown file extension for: %s", file_weights.c_str()));

    if (!file_lengths.empty())
        tl = tvb::csv_load(file_lengths);
    else
        tl = tvb::TArray2d::Zero(C.rows(), C.cols());

    if (file_gaba.ends_with(".csv"))
        gaba_vector = tvb::csv_load(file_weights);
    else if (file_gaba.ends_with(".npy")) {
        gaba_vector = tvb::npy2Vector(file_gaba);
    } else
        throw std::runtime_error(string_format("Unknown file extension for: %s", file_gaba.c_str()));
}


RunParams run(RunParams rp) {

    string f_prefix = (path(rp.path_out) / path(rp.file_prefix)).lexically_normal().string();
    for (auto const &p: rp.params) {
        f_prefix += string_format("_%s_%.2f", p.name.c_str(), p.value);
    }
    string filename = f_prefix + ".png";

    if (!rp.force_output && std::filesystem::exists(filename)) {
        std::cout << string_format("File %s already exists", filename.c_str()) << std::endl;
        delete rp.monitor;
        rp.monitor = nullptr;
        return rp;
    }

    tvb::TArray2d C;
    tvb::TArray2d tl;
    tvb::TArray1d gaba_vector;

    load_data(rp.file_weights, rp.file_lengths, rp.gaba_vector, C, tl, gaba_vector);

    int N = C.rows();

    Float k = 0.15 / (C.rowwise().sum().sum() / N);
    // C *= k;
    // tvb::csv_save("sc_d_norm.csv", C);

    auto *con = new tvb::Connectivity(C, tl, rp.speed);

    milliseconds total_time(0);
    std::cout << string_format("Starting computation for: %s", filename.c_str()) << std::endl;

    // auto *model = new tvb::Montbrio(N, rp.t_start, rp.t_end, rp.dt);
    // auto *model = new tvb::ReducedWongWangExcInh(N);
    auto *model = new ZerlautGABA(N);
    model->configure();
    model->set_param("gaba_ratio", gaba_vector);
    float G = 1.0;
    auto g_it = std::find(rp.params.begin(), rp.params.end(), "G");
    if (g_it != rp.params.end()) {
        G = g_it->value;
        rp.params.erase(g_it);
    }

    for (auto const &p: rp.params)
        if (std::isalpha(p.name[0])) model->set_param(p.name, p.value);

    model->init_dependant();
    // rp.monitor = new tvb::BoldTVB(N, 720.0, rp.dt, {0});
    // rp.monitor = new tvb::BoldBalloonWindkessel(N, 1.0, 720.0, rp.dt, {0});
    // rp.monitor = new tvb::RawSubSample(1.0, rp.dt, {3});

    TArray1d sigmas = TArray1d::Constant(model->n_vars(), 0.0);
    for (auto const &p: rp.params)
        if (p.name[0] == '_') {
            auto idx = std::stoi(p.name.substr(2, 1));
            sigmas[idx] = p.value;
        }


    // sigmas << 3e-5, 3e-5, 0.0, 0.0;
    // auto *integrator = new tvb::EulerStochastic(rp.dt, new Additive(sigmas, rp.dt));
    auto *integrator = new tvb::EulerDeterministic(rp.dt);

    auto coupling = new tvb::CouplingLinearSparse(con->weights(), con->delays(), model->cvars());
    coupling->setScale(G);

    auto start = std::chrono::high_resolution_clock::now();

    SimConfig sim_config;

    sim_config.setModel(model);
    sim_config.setConnectivity(con);
    sim_config.setIntegrator(integrator);
    // sim_config.setMonitor(rp.monitor);
    sim_config.setCoupling(coupling);
    sim_config.setIntegrationInterval(rp.t_start, rp.t_end);
    sim_config.setNumIterations(1);
    sim_config.setDeltaIntegration(0.00001);

    auto [converged, monitor] = tvb::simulate(sim_config, 1.0, 0);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start);

    std::cout << string_format("Computation time (%s): %d msecs", filename.c_str(), duration.count()) << std::endl;

    total_time += duration;

//    std::vector<long unsigned> shape{J_i.rows(), J_i.cols()};
//    npy::SaveArrayAsNumpy(f_prefix + ".npz", false, shape.size(), shape.data(), (Float*)J_i.data());

    // save_fig(monitor, f_prefix);
    save_cvs(monitor, f_prefix);

    delete monitor;
    rp.file_out = filename;

    delete model;
    delete coupling;

    return rp;
}


void to_cout(const std::vector<std::string> &v) {
    std::copy(v.begin(), v.end(), std::ostream_iterator<std::string>{
            std::cout, "\n"});
}

int main(int argc, char **argv) {

    try {
        options_description desc{"Options"};
        desc.add_options()
                ("help,h", "Help screen")
                ("params", value<std::vector<std::string>>()->multitoken()->required(), "Parameters to sweep")
                ("sc-matrix", value<std::string>()->required(), "Structural connectivity matrix")
                ("gaba-vector", value<std::string>()->required(), "Vector with neuroreceptor density")
                ("length-matrix", value<std::string>(), "Connection lengths matrix matrix")
                ("speed", value<float>()->default_value(1e6), "Signal speed")
                ("use-threads", value<bool>()->default_value(false), "Use threads")
                ("time-start", value<float>()->default_value(0.0), "Start of simulation (ms)")
                ("time-end", value<float>()->default_value(10000.0), "End of simulation (ms)")
                ("dt", value<float>()->default_value(0.1), "Integration step (ms)")
                ("force-output", bool_switch()->default_value(false), "Force overwrite output files")
                ("sigmas", value<std::vector<std::string>>()->multitoken(), "Noise sigmas")
                ("var-of-interest", value<int>()->required(), "Variable of interest in the model to explore")
                ("jube-cpu-pp", value<int>()->required(), "Number of cores per execution")
                ("value-base", value<float>()->default_value(base_value),
                 "Point of equilibrium for excitatory intensity output")
                ("out-path", value<std::string>()->required(), "Output path")
                ("algo", value<std::string>()->required(), "Algorithm to run")
                ("norm", value<std::string>()->required(), "SC matrix normalization method")
                ("experiment-name", value<std::string>()->required(), "name for experiment (folder)");
        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);
        notify(vm);

        std::vector<SweepParam> params;
        if (vm.count("help")) {
            std::cout << desc << '\n';
            return 0;
        }
        if (vm.count("params")) {
            for (auto &s: vm["params"].as<std::vector<std::string>>()) {
                if (std::isalpha(s[0])) {
                    if (!params.empty() && !(params.back().values.size() == 1 || params.back().values.size() == 3))
                        throw std::runtime_error(string_format("Malformed parameter <%s>\n", s.c_str()));
                    params.emplace_back();
                    params.back().name = s;
                } else {
                    try {
                        float value = std::stof(s);
                        params.back().values.push_back(value);
                    } catch (const error &ex) {
                        throw std::runtime_error("Syntax error in sweep parameters");
                    }
                }
            }
            if (params.back().values.size() != 1 && params.back().values.size() != 3)
                throw std::runtime_error(string_format("Malformed parameter <%s>\n", params.back().name.c_str()));
        }

        if (vm.count("sigmas")) {
            auto p_size = params.size();
            for (auto &s: vm["sigmas"].as<std::vector<std::string>>()) {
                if (s[0] == '_') {
                    if (!params.empty() && !(params.back().values.size() == 1 || params.back().values.size() == 3))
                        throw std::runtime_error(string_format("Malformed parameter <%s>\n", s.c_str()));
                    params.emplace_back();
                    params.back().name = s;
                } else {
                    try {
                        float value = std::stof(s);
                        params.back().values.push_back(value);
                    } catch (const error &ex) {
                        throw std::runtime_error("Syntax error in sweep parameters");
                    }
                }
            }
            if (params.back().values.size() != 1 && params.back().values.size() != 3)
                throw std::runtime_error(string_format("Malformed parameter <%s>\n", params.back().name.c_str()));

            if (params.size() - p_size != 4)
                throw std::runtime_error(string_format("Wrong number of sigma parameters, it should be 4, its %i",
                                                       params.size() - p_size));
        }


        std::vector<RunParams> param_combs(1);
        for (auto const &p: params) {
            if (p.values.size() == 1) {
                for (auto &pc: param_combs)
                    pc.params.emplace_back(p.name, p.values[0]);
            } else {
                std::vector<RunParams> new_param_combs;
                for (auto v: tvb::range(p.values[0], p.values[1], p.values[2])) {
                    for (auto const &pc: param_combs) {
                        new_param_combs.push_back(pc);
                        new_param_combs.back().params.emplace_back(p.name, v);
                    }
                }
                param_combs = new_param_combs;
            }
        }

        int num_cores = 8;
        if (vm.count("jube-cpu-pp"))
            num_cores = vm["jube-cpu-pp"].as<int>();

        if (param_combs.size() < num_cores) num_cores = param_combs.size();

        if (!vm["use-threads"].as<bool>()) {
            for (auto &pc: param_combs) {
                pc.init(vm);
                run(pc);
            }
        } else {
            tvb::ThreadPool<RunParams> tp(num_cores);
            tp.start();
            for (auto &pc: param_combs) {
                pc.init(vm);
                tp.queue_job([pc] { return run(pc); });
            }

            while (!tp.finished()) {
                std::optional<RunParams> op = tp.get_result();
                sleep(1);
            }

            tp.stop();
        }

    }
    catch (const error &ex) {
        std::cerr << ex.what() << '\n';
    }
}
