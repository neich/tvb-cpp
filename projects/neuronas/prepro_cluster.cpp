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
#include <tvb-cpp/simulator/models/naskar.h>
#include <tvb-cpp/simulator/integrators/euler_stochastic.h>
#include <tvb-cpp/tools/threadpool.h>
#include <tvb-cpp/tools/csv_tools.h>
#include <tvb-cpp/tools/npy.h>
#include <tvb-cpp/simulator/integrators/euler_deterministic.h>
#include <tvb-cpp/tools/observers/sw_fcd.h>
#include <tvb-cpp/tools/observers/ph_fcd.h>
#include <tvb-cpp/tools/bold_filters.h>
#include <tvb-cpp/tools/json.h>
#include <tvb-cpp/simulator/bold/bold_tvb.h>

#include <chrono>
#include <filesystem>
#include <thread>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/process.hpp>
#include <iostream>

#include <utility>
#include <tvb-cpp/tools/algo/load_or_compute.h>

using namespace boost::program_options;
using namespace boost::process;
using namespace std::filesystem;
using namespace tvb;
using namespace std::chrono;
using namespace std;
using json = nlohmann::json;

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
    string file_ts;
    vector<string> norm;
    string model;
    float speed = 1e6;
    float G = 15.0;
    float tr;

    string error;

    RunParams() = default;

    RunParams(std::vector<Parameter> params, tvb::Monitor *monitor) : params(std::move(params)), monitor(monitor) {}

    void init(const variables_map &vm) {
        this->model = vm["model"].as<string>();
        this->voi = vm["var-of-interest"].as<int>();
        this->value_base = vm["value-base"].as<float>();
        this->file_weights = vm["sc-matrix"].as<std::string>();
        if (vm.count("norm"))
            this->norm = vm["norm"].as<vector<string>>();
        this->t_start = vm["time-start"].as<float>();
        this->tr = vm["tr"].as<float>();
        this->t_end = vm["time-end"].as<float>();
        this->dt = vm["dt"].as<float>();
        if (vm.count("length-matrix"))
            this->file_lengths = vm["length-matrix"].as<std::string>();
        if (vm.count("time-series"))
            this->file_ts = vm["time-series"].as<std::string>();
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

string getPrefix(const vector<Parameter> &params);

void collect_results(const string &job_id, const path &out_dir, const vector<RunParams> &param_combs);

void saveJSON(const RunParams &rp, const string &f_prefix, const path &out_dir, double fitting);

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

    tvb::csv_save(filename, y_plot, "", true);

}


TArray2d_uptr load_data(const string &file_weights) {
    if (file_weights.ends_with(".csv"))
        return tvb::csv_load(file_weights);
    else if (file_weights.ends_with(".npz"))
        return tvb::npz2Matrixd(file_weights, "SC");
    else if (file_weights.ends_with(".npy"))
        return tvb::npy2Matrixd(file_weights);
    else
        throw std::runtime_error(string_format("Unknown file extension for: %s", file_weights.c_str()));
}


//RunParams run(RunParams rp, unsigned n = 1, unsigned total = 1) {
//
//
//    string f_prefix = getPrefix(rp.params);
//
//    path out_dir = rp.path_out;
//
//    tvb::TArray2d C;
//    tvb::TArray2d tl;
//
//    load_data(rp.file_weights, rp.file_lengths, C, tl);
//
//    if (rp.norm.size() > 0) {
//        if (rp.norm[0] == "Gus") {
//            auto factor = boost::lexical_cast<double>(rp.norm[1]);
//            double maxC = C.rowwise().sum().maxCoeff();
//            C = C / maxC * factor;
//        } else {
//            throw std::runtime_error(string_format("Unsupported normalization method <%s>\n", rp.norm[0].c_str()));
//        }
//    }
//
//    int N = C.rows();
//
//    Float k = 0.15 / (C.rowwise().sum().sum() / N);
//    // C *= k;
//    // tvb::csv_save("sc_d_norm.csv", C);
//
//    auto *con = new tvb::Connectivity(C, tl, rp.speed);
//
//    milliseconds total_time(0);
//    std::cout << string_format("Starting computation (%d of %d)for: %s", n, total, f_prefix.c_str()) << std::endl;
//
//    // auto *model = new tvb::Montbrio(N, rp.t_start, rp.t_end, dt);
//    // auto *model = new tvb::ReducedWongWangExcInh(N);
//    tvb::Model *model;
//    TArray1d sigmas;
//    if (rp.model == "ZerlautGABA") {
//        model = new ZerlautGABA(N);
//        sigmas = TArray1d::Constant(model->n_vars(), 0.0);
//        sigmas[model->n_vars() - 1] = 1.0;
//        model->configure();
//        if (gaba_vector.size() > 0)
//            model->set_param("gaba_ratio", gaba_vector);
//        else
//            model->set_param("gaba_ratio", TArray1d::Ones(N));
//    } else if (rp.model == "Montbrio") {
//        model = new Montbrio(N);
//        sigmas = TArray1d::Constant(model->n_vars(), 0.0);
//        model->configure();
//    }
//    else {
//        throw std::runtime_error(string_format("Unknown model <%s>\n", rp.model.c_str()));
//    }
//
//    float G = 1.0;
//    auto g_it = std::find(rp.params.begin(), rp.params.end(), "G");
//    if (g_it != rp.params.end()) {
//        G = g_it->value;
//        // rp.params.erase(g_it);
//    }
//
//    for (auto const &p: rp.params)
//        if (std::isalpha(p.name[0]) && p.name != "G") model->set_param(p.name, p.value);
//
//    model->init_dependant();
//    // rp.monitor = new tvb::BoldTVB(N, 720.0, dt, {0});
//    // rp.monitor = new tvb::BoldBalloonWindkessel(N, 1.0, 720.0, dt, {0});
//    // rp.monitor = new tvb::RawSubSample(1.0, dt, {3});
//
///*
//    for (auto const &p: rp.params)
//        if (p.name[0] == '_') {
//            auto idx = std::stoi(p.name.substr(2, 1));
//            sigmas[idx] = p.value;
//        }
//*/
//
//
//    sigmas << 0, 0, 0, 0, 0, 0, 0, 1;
//    auto *integrator = new tvb::EulerStochastic(dt, new Additive(sigmas, dt));
//    // auto *integrator = new tvb::EulerDeterministic(dt);
//
//    auto coupling = new tvb::CouplingLinearSparse(con->weights(), con->delays(), model->cvars());
//    coupling->setScale(G);
//
//    if (rp.algo == "explore_G") {
//        path npy_file = out_dir;
//        npy_file /= rp.job_id + f_prefix + ".npy";
//
//        if (!rp.force_output && std::filesystem::exists(npy_file)) {
//            std::cout << string_format("File %s already exists", npy_file.c_str()) << std::endl;
//            delete rp.monitor;
//            rp.monitor = nullptr;
//            return rp;
//        }
//
//        auto start = std::chrono::high_resolution_clock::now();
//
//        SimConfig sim_config;
//
//        sim_config.setModel(model);
//        sim_config.setConnectivity(con);
//        sim_config.setIntegrator(integrator);
//        // sim_config.setMonitor(rp.monitor);
//        sim_config.setCoupling(coupling);
//        sim_config.setIntegrationInterval(rp.t_start, rp.t_end);
//        sim_config.setNumIterations(1);
//        sim_config.setDeltaIntegration(0.00001);
//
//        Simulator simulator{};
//        TArray2d initial_state = TArray2d::Zero(C.cols(), model->n_vars());
//        auto *monitor = new TemporalAverage(con->weights().cols(), rp.ta_period, dt, {0});
//
//        simulator.run(sim_config.model(),
//                      sim_config.connectivity(),
//                      sim_config.integrator(),
//                      {monitor},
//                      sim_config.coupling(),
//                      0, rp.t_end,
//                      nullptr,
//                      &initial_state);
//
//        auto stop = std::chrono::high_resolution_clock::now();
//        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
//                std::chrono::high_resolution_clock::now() - start);
//
//        // std::cout << string_format("Computation time (%s): %d msecs", npy_file.c_str(), duration.count()) << std::endl;
//        std::cout << string_format("Computation time: %d msecs", duration.count()) << std::endl;
//
//        total_time += duration;
//
//        // save_fig(monitor, f_prefix);
//        TArray2d data = monitor->voi2Array(0);
//        // Matrixd2np(data, npy_file.string());
//
//        delete monitor;
//        // rp.file_out = npy_file.string();
//
//        delete model;
//        delete coupling;
//    } else {
//
//        path npz_file = out_dir;
//        // npz_file /= rp.job_id + f_prefix + ".npz";
//        npz_file /= f_prefix + ".npz";
//
//        if (std::filesystem::exists(npz_file) && !rp.force_output) {
//            std::cout << string_format("File %s already exists\n", npz_file.c_str()) << std::flush;
//            auto data = npz2MatrixdMap(npz_file.string());
//            saveJSON(rp, f_prefix, out_dir, data["fit"](0,0));
//            delete rp.monitor;
//            rp.monitor = nullptr;
//            delete model;
//            delete coupling;
//            return rp;
//        }
//
//        path pe_file = out_dir / "fNeuro_emp.npy";
//        if (!exists(pe_file))
//            throw std::runtime_error(string_format("Preprocessed file does not exists %s\n", pe_file.c_str()));
//
//        TArray2dMap data = npz2MatrixdMap(pe_file.string());
//        TArray2d processed_emp = data["swFCD"];
//
//        BandPassFilter bpf(0.008, 0.08, 2.5);
//        // SW_FC measure(30, 10, true, bpf);
//        PhFCD measure(5, true, bpf);
//
//        int N = data["nsub"](0, 0);
//        measure.init(N, N);
//        auto start = std::chrono::high_resolution_clock::now();
//
//        SimConfig sim_config;
//
//        sim_config.setModel(model);
//        sim_config.setConnectivity(con);
//        sim_config.setIntegrator(integrator);
//        sim_config.setCoupling(coupling);
//        sim_config.setNumIterations(1);
//        sim_config.setDeltaIntegration(0.00001);
//
//        Simulator simulator{};
//        TemporalAverage* ta_mon = new TemporalAverage(N, 1, dt, {0});
//        TArray2d initial_state = TArray2d::Zero(C.cols(), model->n_vars());
//
//        simulator.run(sim_config.model(),
//                      sim_config.connectivity(),
//                      sim_config.integrator(),
//                      {ta_mon},
//                      sim_config.coupling(),
//                      0, rp.t_end,
//                      nullptr,
//                      &initial_state);
//
//
//        TArray2d raw_signal = ta_mon->voi2Array(0);
//        BoldTVB *btvb = new BoldTVB(TR);
//        auto [bold_times, bold_signal] = btvb->compute_bold(raw_signal, 1.0);
//        TArray2d proc_signal = measure.from_fMRI(bold_signal);
//        measure.accumulate(proc_signal);
//
//
//        auto measureValues = measure.postprocess();
//        auto fitting = measure.distance(measureValues, processed_emp);
//        cout << string_format("Distance for <%s> : <%f>\n", f_prefix.c_str(), fitting);
//
//        TArray2dMap npz_data;
//        npz_data["measure"] = measureValues;
//        npz_data["fit"] = {{(double) fitting}};
//        MatrixdMap2npz(npz_file.string(), npz_data);
//
//        auto stop = std::chrono::high_resolution_clock::now();
//        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
//                std::chrono::high_resolution_clock::now() - start);
//        cout << string_format("Computed time series (%d, %d) for <%s> (time: <%d>)\n", bold_signal.rows(),
//                              bold_signal.cols(), f_prefix.c_str(), duration.count()) << flush;
//
//        saveJSON(rp, f_prefix, out_dir, fitting);
//    }
//
//    return rp;
//}

void saveJSON(const RunParams &rp, const string &f_prefix, const path &out_dir, double fitting) {
//    auto json_file = out_dir / (rp.job_id + "_" + f_prefix + ".json");
//    ofstream jsonf(json_file, std::ios_base::out | std::ios_base::trunc);
//    json oj;
//    oj["fit"] = fitting;
//    json pj;
//    for (auto const &p: rp.params) {
//        pj[p.name.c_str()] = p.value;
//    }
//    oj["params"] = pj;
//    jsonf << oj;
//    jsonf.close();
}

string getPrefix(const vector<Parameter> &params) {
    string f_prefix;
    for (auto const &p: params) {
        if (f_prefix.size() > 0)
            f_prefix += "_";
        f_prefix += string_format("%s_%.2f", p.name.c_str(), p.value);
    }
    return f_prefix;
}

TArray2d simulateSingleSubject(const RunParams &params, SW_FC &fc) {
    return TArray2d{};
}


void to_cout(const std::vector<std::string> &v) {
    std::copy(v.begin(), v.end(), std::ostream_iterator<std::string>{
            std::cout, "\n"});
}

TArray2d_uptr processBOLDSignals(const tvb::TArray2dMap &bolds) {
    int num_subjects = bolds.size();
    int n_roi = bolds.begin()->second->rows();

    BandPassFilter bpf(0.01, 0.1, 2.0);
    PhFCD measure(10, true, bpf);

    measure.init(num_subjects, n_roi);

    auto it = bolds.begin();
    for (int pos = 0; pos < num_subjects; ++pos, ++it) {
        printff("Processing BOLD subject %i\n", pos);
        TArray2d procSignal = measure.from_fMRI(*it->second);
        measure.accumulate(procSignal, pos);
    }
    return std::move(measure.postprocess());
}

int main(int argc, char **argv) {

        options_description desc{"Options"};
        desc.add_options()
                ("help,h", "Help screen")
                ("tr", value<float>()->default_value(2000.0), "Bold sampling rate (ms)")
                ("we", value<float>(), "G value to compute")
                ("we-range", value<std::vector<float>>()->multitoken(), "G range to axplore")
                ("force-output", bool_switch()->default_value(false), "Force overwriting of CSV files if they already exists")
                ("out-path", value<std::string>()->required(), "Output path")
                ("data-path", value<std::string>(), "Data path");

        variables_map vm;
        store(command_line_parser(argc, argv)
                      .options(desc)
                      .style(command_line_style::unix_style ^ command_line_style::allow_short)
                      .run(), vm);
        notify(vm);

        if (vm.count("help")) {
            std::cout << desc << '\n';
            return 0;
        }

        if (vm.count("we")) {
            path out_path = vm["out-path"].as<string>();
            path data_path = vm["data-path"].as<string>();

            if (!exists(out_path)) {
                create_directory(out_path);
                printff("Output directory created: %s\n", out_path.c_str());
            }

            auto sc_filename = data_path / "SC_dbs80HARDIFULL.npz";

            auto sc_data = tvb::npz2MatrixdMap(sc_filename);
            tvb::TArray2d_uptr C = tvb::npz2Matrixd(sc_filename, "SC");
            tvb::TArray2d tl = tvb::TArray2d::Ones(C->rows(), C->cols());

            int N = C->rows();

            Float k = 1.0 / (C->maxCoeff() * 0.1);
            (*C) *= k;

            auto *con = new tvb::Connectivity(*C, tl, 1e6); // TODO: Connectivity should probably use shared_ptr

            milliseconds total_time(0);
            // std::cout << string_format("Starting computation (%d of %d)for: %s", n, total, f_prefix.c_str()) << std::endl;

            // auto *model = new tvb::Montbrio(N, rp.t_start, rp.t_end, dt);
            // auto *model = new tvb::ReducedWongWangExcInh(N);
            tvb::Model *model = new tvb::Naskar(N);
            model->configure();
            model->init_dependant();

            TArray1d sigmas(3);
            sigmas << 1e-4, 0, 0;
            tvb::Float dt = 0.1;
            auto *integrator = new tvb::EulerStochastic(dt, new Additive(sigmas, dt));

            auto coupling = new tvb::CouplingLinearSparse(con->weights(), con->delays(), model->cvars());
            auto we = vm["we"].as<float>();
            coupling->setScale(we);



            BandPassFilter bpf(0.01, 0.1, 2.0);
            PhFCD measure(10, true, bpf);

            measure.init(N, N);
            cout << "Processing empirical subjects ..." << endl;
            auto proc_emp_file = out_path / "fNeuroEmp_REST.npz";

            TArray2d_uptr processed_emp;
            if (exists(proc_emp_file)) {
                cout << "Loading file: " << proc_emp_file << endl;
                processed_emp = npz2Matrixd(proc_emp_file, "phFCD");
            } else {
                auto kk = string_format("Computing file: %s\n", proc_emp_file.c_str());
                printff("Computing file: %s\n", proc_emp_file.c_str());
                path fmris_file = data_path / "hcp1003_REST_LR_dbs80.npz";
                auto fmris = tvb::npz2MatrixdMap(fmris_file);
                processed_emp = processBOLDSignals(fmris);
                cout << "Finished processing subjects!" << endl << flush;
                tvb::MatrixdMap2npz(proc_emp_file, {{"phFCD", std::move(processed_emp)}});
            }

            auto start = std::chrono::high_resolution_clock::now();

            SimConfig sim_config;

            sim_config.setModel(model);
            sim_config.setConnectivity(con);
            sim_config.setIntegrator(integrator);
            sim_config.setCoupling(coupling);
            sim_config.setNumIterations(1);
            sim_config.setDeltaIntegration(0.00001);

            Simulator simulator{};
            TemporalAverage *ta_mon = new TemporalAverage(N, 1, dt, {0});
            TArray2d initial_state = TArray2d::Zero(C->cols(), model->n_vars());

            auto TR = 2.0;
            auto dtt = 1e-3;
            auto tmax = 220.0*TR/dtt;
            printff("Starting simulation for we = %f.2\n", we);
            simulator.run(sim_config.model(),
                          sim_config.connectivity(),
                          sim_config.integrator(),
                          {ta_mon},
                          sim_config.coupling(),
                          0, tmax,
                          nullptr,
                          &initial_state);


            TArray2d raw_signal = ta_mon->voi2Array(0);
            BoldTVB *btvb = new BoldTVB(TR*1000.0);
            printff("Computing BOLD for we = %f.2\n", we);
            auto [bold_times, bold_signal] = btvb->compute_bold(raw_signal, 1.0);
            printff("From fMRI ... \n");
            TArray2d proc_signal = measure.from_fMRI(bold_signal.transpose());
            printff("Accumulate ... \n");
            measure.accumulate(proc_signal);

            printff("Postprocess ... \n");
            auto measureValues = measure.postprocess();
            printff("Distance ... \n");
            auto fitting = measure.distance(*measureValues, *processed_emp);
            printff("Distance for we = %f.2 : %f\n", we, fitting);

            TArray2dMap npz_data;
            auto npz_file = out_path / string_format("fitting_we_%f.2.npz", we);
            npz_data["measure"] = std::move(measureValues);
            npz_data["fit"] = std::make_shared<TArray2d>(TArray2d({{(tvb::Float) fitting}}));
            MatrixdMap2npz(npz_file.string(), npz_data);

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::high_resolution_clock::now() - start);
            printff("Computed time series (%d, %d) for we = <%f.2> (time: <%d>)\n", bold_signal.rows(),
                                  bold_signal.cols(), we, duration.count());

            // saveJSON(rp, f_prefix, out_dir, fitting);
        }
        if (vm["srun"].as<bool>()) {
            cout << "Running jobs using srun/slurm" << endl;
            std::vector<child *> processes;
            auto srun = search_path("srun");
            unsigned srun_pack_size = 100;
                string args = "";
                args += string_format(" -N 1 -n 1 -c 1 --exclusive --mem-per-cpu=2000MB %s",
                                      std::filesystem::canonical("/proc/self/exe").c_str());
                args += string_format(" --sc-matrix %s", vm["sc-matrix"].as<string>().c_str());
                args += string_format(" --out-path %s", vm["out-path"].as<string>().c_str());
                args += string_format(" --algo %s", vm["algo"].as<string>().c_str());
                args += string_format(" --experiment-name %s", vm["experiment-name"].as<string>().c_str());
                args += string_format(" --job-id %s", vm["job-id"].as<string>().c_str());
                if (vm.count("norm")) {
                    args += " --norm";
                    for (auto &s: vm["norm"].as<vector<string>>())
                        args += " " + s;
                }
                if (vm.count("gaba-vector") > 0)
                    args += string_format(" --gaba-vector %s", vm["gaba-vector"].as<string>().c_str());
                if (vm["force-output"].as<bool>())
                    args += " --force-output";
                if (vm.count("tr") > 0)
                    args += string_format(" --tr %f", vm["tr"].as<float>());
                if (vm.count("model") > 0)
                    args += string_format(" --model %s", vm["model"].as<string>().c_str());
                if (vm.count("time-series") > 0)
                    args += string_format(" --time-series %s", vm["time-series"].as<string>().c_str());
                if (vm.count("ta-period") > 0)
                    args += string_format(" --ta-period %f", vm["ta-period"].as<float>());
                if (vm.count("time-end") > 0)
                    args += string_format(" --time-end %f", vm["time-end"].as<float>());
                if (vm.count("time-start") > 0)
                    args += string_format(" --time-start %f", vm["time-start"].as<float>());


                cout << string_format("Executing [%s, %s]\n", srun.c_str(), args.c_str()) << std::endl;
                auto *c = new child(srun.string() + args);
                processes.push_back(c);
                if (processes.size() == srun_pack_size) {
                    // cout << string_format("Waiting for pack of srun processes to finish (%d total sent)\n", n);
                    for (auto c: processes)
                        c->wait();
                    processes.clear();
                }


            for (auto c: processes)
                c->wait();

            // cout << string_format("Collecting results for %d simulations\n", param_combs.size());
            // collect_results(vm["job-id"].as<string>(), out_path, param_combs);

        }
}

void collect_results(const string &job_id, const path &out_dir, const vector<RunParams> &param_combs) {
    auto json_file = out_dir / (job_id + ".json");
    ofstream jsonf(json_file, std::ios_base::out | std::ios_base::trunc);
    json oj;
    for (auto &pc: param_combs) {
        string f_prefix = getPrefix(pc.params);
        try {
            auto sim_file = out_dir / (job_id + "_" + f_prefix + ".json");
            ifstream comb_file(sim_file);
            json cjson;
            comb_file >> cjson;
            oj[f_prefix] = cjson;
        } catch (const std::exception& e) {
            std::cout << "Error parsing JSON file: " << json_file << std::endl;
            std::cout << "Reason: " << e.what() << std::endl;
            return;
        }
    }
    jsonf << oj;
    jsonf.close();
}
