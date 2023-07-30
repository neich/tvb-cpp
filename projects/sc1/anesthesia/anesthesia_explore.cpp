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
#include <tvb-cpp/simulator/simulator.h>
#include <tvb-cpp/simulator/monitor.h>
#include <tvb-cpp/simulator/models/reduced_ww_ext.h>
#include <tvb-cpp/simulator/models/montbrio.h>
#include <tvb-cpp/simulator/integrators/euler_stochastic.h>
#include <tvb-cpp/tools/threadpool.h>
#include <tvb-cpp/tools/csv_tools.h>
#include <tvb-cpp/tools/npy.h>
#include <tvb-cpp/simulator/integrators/euler_deterministic.h>
#include <tvb-cpp/tools/algo/fic/functions/observers/sw_fcd.h>
#include <tvb-cpp/tools/algo/fic/functions/bold_filters.h>
#include <tvb-cpp/simulator/monitors/bold_tvb.h>
#include "zerlaut_gaba.h"

#include <chrono>
#include <filesystem>
#include <thread>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/process.hpp>
#include <iostream>

#include <utility>

using namespace boost::program_options;
using namespace boost::process;
using namespace std::filesystem;
using namespace tvb;
using namespace std::chrono;
using namespace std;

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
    string experiment_name;
    string file_prefix;
    string path_out;
    string file_weights;
    string file_lengths;
    string file_ts;
    string gaba_vector;
    string algo;
    vector<string> norm;
    string job_id;
    float speed = 1e6;
    float G = 1.0;

    RunParams() = default;

    RunParams(std::vector<Parameter> params, tvb::Monitor *monitor) : params(std::move(params)), monitor(monitor) {}

    void init(const variables_map &vm) {
        this->algo = vm["algo"].as<string>();
        this->voi = vm["var-of-interest"].as<int>();
        this->value_base = vm["value-base"].as<float>();
        this->file_weights = vm["sc-matrix"].as<std::string>();
        if (vm.count("norm"))
            this->norm = vm["norm"].as<vector<string>>();
        if (vm.count("gaba-vector"))
            this->gaba_vector = vm["gaba-vector"].as<std::string>();
        this->t_start = vm["time-start"].as<float>();
        this->t_end = vm["time-end"].as<float>();
        this->dt = vm["dt"].as<float>();
        if (vm.count("length-matrix"))
            this->file_lengths = vm["length-matrix"].as<std::string>();
        if (vm.count("time-series"))
            this->file_ts = vm["time-series"].as<std::string>();
        this->job_id = vm["job-id"].as<std::string>();
        // this->file_prefix = vm["out-file-prefix"].as<std::string>();
        this->path_out = vm["out-path"].as<std::string>();
        this->experiment_name = vm["experiment-name"].as<std::string>();
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

void save_cvs(tvb::Monitor *monitor, const string &filename) {
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

    tvb::csv_save(filename, y_plot, true);

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

    if (file_gaba.size() == 0)
        gaba_vector = TArray1d::Ones(C.cols());
    else {
        if (file_gaba.ends_with(".csv"))
            gaba_vector = tvb::csv_load(file_weights);
        else if (file_gaba.ends_with(".npy")) {
            gaba_vector = tvb::npy2Vector(file_gaba);
        } else
            throw std::runtime_error(string_format("Unknown file extension for: %s", file_gaba.c_str()));
    }
}


RunParams run(RunParams rp) {

    string f_prefix;
    for (auto const &p: rp.params) {
        f_prefix += string_format("_%s_%.2f", p.name.c_str(), p.value);
    }

    path out_dir = rp.path_out;
    out_dir /= rp.experiment_name;

    tvb::TArray2d C;
    tvb::TArray2d tl;
    tvb::TArray1d gaba_vector;

    load_data(rp.file_weights, rp.file_lengths, rp.gaba_vector, C, tl, gaba_vector);

    if (rp.norm.size() > 0) {
        if (rp.norm[0] == "Gus") {
            auto factor = boost::lexical_cast<double>(rp.norm[1]);
            double maxC = C.rowwise().sum().maxCoeff();
            C = C / maxC * factor;
        } else {
            throw std::runtime_error(string_format("Unsupported normalization method <%s>\n", rp.norm[0].c_str()));
        }
    }

    int N = C.rows();

    Float k = 0.15 / (C.rowwise().sum().sum() / N);
    // C *= k;
    // tvb::csv_save("sc_d_norm.csv", C);

    auto *con = new tvb::Connectivity(C, tl, rp.speed);

    milliseconds total_time(0);
    std::cout << string_format("Starting computation for: %s", f_prefix.c_str()) << std::endl;

    // auto *model = new tvb::Montbrio(N, rp.t_start, rp.t_end, rp.dt);
    // auto *model = new tvb::ReducedWongWangExcInh(N);
    auto *model = new ZerlautGABA(N);
    model->configure();
    if (gaba_vector.size() > 0)
        model->set_param("gaba_ratio", gaba_vector);
    else
        model->set_param("gaba_ratio", TArray1d::Ones(N));

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
    sigmas[model->n_vars()-1] = 1.0;
/*
    for (auto const &p: rp.params)
        if (p.name[0] == '_') {
            auto idx = std::stoi(p.name.substr(2, 1));
            sigmas[idx] = p.value;
        }
*/


    // sigmas << 3e-5, 3e-5, 0.0, 0.0;
    auto *integrator = new tvb::EulerStochastic(rp.dt, new Additive(sigmas, rp.dt));
    // auto *integrator = new tvb::EulerDeterministic(rp.dt);

    auto coupling = new tvb::CouplingLinearSparse(con->weights(), con->delays(), model->cvars());
    coupling->setScale(G);

    if (rp.algo == "explore_G") {
        path cvs_file = out_dir;
        cvs_file /= rp.job_id + f_prefix + ".csv";

        if (!rp.force_output && std::filesystem::exists(cvs_file)) {
            std::cout << string_format("File %s already exists", cvs_file.c_str()) << std::endl;
            delete rp.monitor;
            rp.monitor = nullptr;
            return rp;
        }

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

        Simulator simulator{};
        TArray2d initial_state = TArray2d::Zero(C.cols(), model->n_vars());
        auto *monitor = new TemporalAverage(con->weights().cols(), 1.0, 0.1, {0});

        simulator.run(sim_config.model(),
                      sim_config.connectivity(),
                      sim_config.integrator(),
                      {monitor},
                      sim_config.coupling(),
                      0, 10000,
                      nullptr,
                      &initial_state);

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - start);

        std::cout << string_format("Computation time (%s): %d msecs", cvs_file.c_str(), duration.count()) << std::endl;

        total_time += duration;

        // save_fig(monitor, f_prefix);
        save_cvs(monitor, cvs_file);

        delete monitor;
        rp.file_out = cvs_file;

        delete model;
        delete coupling;
    } else {

        path pe_file = out_dir / "fNeuro_emp.npy";
        if (!exists(pe_file))
            throw std::runtime_error(string_format("Preprocessed file does not exists %s\n", pe_file.c_str()));

        TArray2dMap data = npz2MatrixdMap(pe_file.c_str());
        TArray2d processed_emp = data["swFCD"];

        BandPassFilter bpf(0.008, 0.08, 2.4);
        SW_FC measure(80, 18, true, bpf);

        int N = data["nsub"](0, 0);
        measure.init(N, N);
        for (int n = 0; n < N; ++n) {
            auto start = std::chrono::high_resolution_clock::now();

            SimConfig sim_config;

            sim_config.setModel(model);
            sim_config.setConnectivity(con);
            sim_config.setIntegrator(integrator);
            sim_config.setCoupling(coupling);
            sim_config.setIntegrationInterval(rp.t_start, rp.t_end);
            sim_config.setNumIterations(1);
            sim_config.setDeltaIntegration(0.00001);

            Simulator simulator{};
            BoldTVB *btvb = new BoldTVB(con->weights().cols(), 2500.0, 0.1, {0});
            TArray2d initial_state = TArray2d::Zero(C.cols(), model->n_vars());

            simulator.run(sim_config.model(),
                          sim_config.connectivity(),
                          sim_config.integrator(),
                          {btvb},
                          sim_config.coupling(),
                          0, 48000,
                          nullptr,
                          &initial_state);

            TArray2d bold_signal = btvb->voi2Array(0);
            TArray2d proc_signal = measure.from_fMRI(bold_signal);
            measure.accumulate(proc_signal);


            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::high_resolution_clock::now() - start);

            total_time += duration;
        }
        auto measureValues = measure.postprocess();
        measure.distance(measureValues, processed_emp);
    }

    return rp;
}

TArray2d simulateSingleSubject(const RunParams &params, SW_FC &fc) {
}


void to_cout(const std::vector<std::string> &v) {
    std::copy(v.begin(), v.end(), std::ostream_iterator<std::string>{
            std::cout, "\n"});
}

TArray2d processBOLDSignals(const vector<TArray2d> &bolds, SW_FC &measure) {
    int NumSubjects = bolds.size();
    int N = bolds[0].rows();

    measure.init(NumSubjects, N);

    for (int pos = 0; pos < bolds.size(); ++pos) {
        TArray2d procSignal = measure.from_fMRI(bolds[pos]);
        measure.accumulate(procSignal, pos);
    }
    return measure.postprocess();
}

int main(int argc, char **argv) {

    try {
        options_description desc{"Options"};
        desc.add_options()
                ("help,h", "Help screen")
                ("param", value<std::vector<std::string>>()->multitoken()->required(), "Parameters to sweep")
                ("norm", value<std::vector<std::string>>()->multitoken(), "Matrix normalilzation method")
                ("sc-matrix", value<std::string>()->required(), "Structural connectivity matrix")
                ("gaba-vector", value<std::string>(), "Vector with neuroreceptor density")
                ("length-matrix", value<std::string>(), "Connection lengths matrix matrix")
                ("speed", value<float>()->default_value(1e6), "Signal speed")
                ("use-threads", bool_switch()->default_value(false), "Use threads")
                ("srun", value<bool>()->default_value(false), "Use srun for parallelism")
                ("time-start", value<float>()->default_value(0.0), "Start of simulation (ms)")
                ("time-end", value<float>()->default_value(10000.0), "End of simulation (ms)")
                ("dt", value<float>()->default_value(0.1), "Integration step (ms)")
                ("force-output", bool_switch()->default_value(false), "Force overwrite output files")
                ("sigmas", value<std::vector<std::string>>()->multitoken(), "Noise sigmas")
                ("var-of-interest", value<int>()->default_value(0), "Variable of interest in the model to explore")
                ("jube-cpu-pp", value<int>()->default_value(1), "Number of cores per execution")
                ("value-base", value<float>()->default_value(base_value),
                 "Point of equilibrium for excitatory intensity output")
                ("out-path", value<std::string>()->required(), "Output path")
                ("time-series", value<std::string>(), "BOLD time series (numSubjects x Time x ROI)")
                ("algo", value<std::string>()->required(), "Algorithm to run")
                ("job-id", value<std::string>()->default_value("default"), "SC matrix normalization method")
                ("experiment-name", value<std::string>()->required(), "name for experiment (folder)");
        variables_map vm;
        store(command_line_parser(argc, argv)
                      .options(desc)
                      .style(command_line_style::unix_style ^ command_line_style::allow_short)
                      .run(), vm);
        notify(vm);

        path out_dir = vm["out-path"].as<string>();
        if (!exists(out_dir))
            throw std::runtime_error(string_format("Output directory does not exists: %s\n", out_dir.c_str()));

        out_dir /= vm["experiment-name"].as<string>();

        if (!exists(out_dir))
            create_directory(out_dir);

        std::vector<SweepParam> params;
        if (vm.count("help")) {
            std::cout << desc << '\n';
            return 0;
        }
        if (vm.count("param")) {
            for (auto &s: vm["param"].as<std::vector<std::string>>()) {
                if (std::isalpha(s[0])) {
                    if (!params.empty() && !(params.back().values.size() == 1 || params.back().values.size() == 3))
                        throw std::runtime_error(string_format("Malformed parameter %s\n", s.c_str()));
                    params.emplace_back();
                    params.back().name = s;
                } else {
                    try {
                        float value = std::stof(s);
                        params.back().values.push_back(value);
                    } catch (...) {
                        throw std::runtime_error(string_format("Syntax error in parameter value: <%s>", s.c_str()));
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
                    } catch (...) {
                        throw std::runtime_error(string_format("Syntax error in sigma value: <%s>", s.c_str()));
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

        if (vm["algo"].as<string>() == "fitting") {
            if (vm.count("time-series") == 0)
                throw std::runtime_error("No time series BOLD provided for algorithm fitting\n");

            path pe_file = vm["out-path"].as<string>();
            pe_file /= vm["experiment-name"].as<string>();
            pe_file /= "fNeuro_emp.npy";
            if (!exists(pe_file)) {

                vector<TArray2d> ts = np2VecMatrixd(vm["time-series"].as<string>());
                BandPassFilter bpf(0.008, 0.08, 2.4);
                SW_FC measure(80, 18, true, bpf);
                vector<TArray2d> transformed_ts;
                for (auto &b: ts)
                    transformed_ts.emplace_back(b.transpose());
                TArray2d processed_emp = processBOLDSignals(transformed_ts, measure);
                TArray2dMap data;
                data["swFCD"] = processed_emp;
                data["nsub"] = (TArray2d(1, 1) << ts.size()).finished();
                data["nsamples"] = (TArray2d(1, 1) << ts[0].rows()).finished();
                MatrixdMap2npz(pe_file.c_str(), data);
            }
        }

        if (vm["srun"].as<bool>()) {
            std::vector<child> processes;
            processes.reserve(param_combs.size());
            auto srun = search_path("srun");
            for (auto &pc: param_combs) {
                string args = "";
                args += string_format(" -N 1 -n 1 -c 1 --exclusive --mem-per-cpu=2000MB %s",
                                      std::filesystem::canonical("/proc/self/exe").c_str());
                args += string_format(" --sc-matrix %s", vm["sc-matrix"].as<string>().c_str());
                args += string_format(" --out-path %s", vm["out-path"].as<string>().c_str());
                args += string_format(" --algo %s", vm["algo"].as<string>().c_str());
                args += string_format(" --experiment-name %s", vm["experiment-name"].as<string>().c_str());
                args += string_format(" --job-id %s", vm["job-id"].as<string>().c_str());
                args += string_format(" --norm %s", vm["norm"].as<string>().c_str());
                if (vm.count("gaba-vector") > 0)
                    args += string_format(" --gaba-vector %s", vm["gaba-vector"].as<string>().c_str());

                for (auto const &p: pc.params)
                    args += string_format(" --param %s %f", p.name.c_str(), p.value);

                processes.emplace_back(srun, args);
            }

            for (auto &c: processes)
                c.wait();

        } else {
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
    }
    catch (const std::runtime_error &ex) {
        std::cerr << ex.what() << '\n';
    }
}
