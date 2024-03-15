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
#include <tvb-cpp/tools/csv_tools.h>
#include <tvb-cpp/tools/npy.h>
#include <tvb-cpp/simulator/integrators/euler_deterministic.h>
#include <tvb-cpp/simulator/models/zerlaut.h>
#include <tvb-cpp/simulator/simulator.h>
#include <tvb-cpp/simulator/monitors/bold_tvb.h>
#include <tvb-cpp/simulator/monitors/bold_BalloonWindkessel.h>
#include <tvb-cpp/tools/algo/fic/functions/balance_fic.h>

#include <tvb-cpp/matplotlibcpp.h>

#include <filesystem>
#include <thread>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <utility>

using namespace boost::program_options;
using namespace std::filesystem;
using namespace tvb;
using namespace std::chrono;
namespace plt = matplotlibcpp;

float base_value = 0.3772258064;
std::mutex matplotlib_mutex;                  // Prevents data races to the job queue

struct SweepParam {
    string name;
    std::vector<float> values;
};

struct Parameter {
    string name;
    float value;

    Parameter(string name, float value) : name(std::move(name)), value(value) {}
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
    float speed = 1e6;

    RunParams() = default;

    RunParams(std::vector<Parameter> params, tvb::Monitor *monitor) : params(std::move(params)), monitor(monitor) {}
};


void save_fig(tvb::Monitor* monitor, const string& file_prefix) {
    if (monitor == nullptr) return;

    int n_records = monitor->getRecords().size();
    int N = monitor->getRecords()[0].record.rows();
    std::vector<std::vector<Float>> y_plot(N, std::vector<Float>(n_records));
    for (unsigned t = 0; t < n_records; ++t) {
        const Monitor::Record &r = monitor->getRecords()[t];
        for (unsigned n = 0; n < N; ++n)
            y_plot[n][t] = r.record(n, 0);
    }

    // tvb::csv_save("./test_simulationRWW_TVB_CPP.csv", y_plot);

    // Plot line from given x and y data. Color is selected automatically.
    std::vector<Float> ls(n_records);
    std::transform(monitor->getRecords().begin(), monitor->getRecords().end(), ls.begin(),
                   [](const Monitor::Record &r) { return r.time/1000; });

    {
        std::unique_lock<std::mutex> lock(matplotlib_mutex);

        for (unsigned n = 0; n < N; ++n) {
            plt::plot(ls, y_plot[n]);
        }
        // Plot a red dashed line from given x and y data.
        // plt::plot(x, w,"r--");
        // Plot a line whose name will show up as "log(x)" in the legend.

        string title = "Reduced Wong Wang";

        plt::title(title);
        plt::ylabel("Ie");
        plt::xlabel("Seconds");
        plt::ylim(0.35, 0.4);
        plt::axhline(base_value - 0.005);
        plt::axhline(base_value + 0.005);
        // Save the image (file format is determined by the extension)
        plt::save(file_prefix + ".png", 300);

        plt::clf();
    }
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
    if (rp.file_weights.ends_with(".csv"))
        C = tvb::csv_load(rp.file_weights);
    else if (rp.file_weights.ends_with(".npz"))
        C = tvb::npz2Matrixd(rp.file_weights, "SC");
    else
        throw std::runtime_error(string_format("Unknown file extension for: %s", rp.file_weights.c_str()));

    int N = C.rows();

    Float k = 0.15 / (C.rowwise().sum().sum() / N);
    // C *= k;
    // tvb::csv_save("sc_d_norm.csv", C);

    tvb::TArray2d tl;
    if (!rp.file_lengths.empty())
        tl = tvb::csv_load(rp.file_lengths);
    else
        tl = tvb::TArray2d::Zero(C.rows(), C.cols());
    auto *con = new tvb::Connectivity(C, tl, rp.speed);

    milliseconds total_time(0);
    std::cout << string_format("Starting computation for: %s", filename.c_str()) << std::endl;

    // auto *model = new tvb::Montbrio(N, rp.t_start, rp.t_end, rp.dt);
    auto *model = new tvb::ReducedWongWangExcInh(N);
    // auto *model = new tvb::ZerlautAdptationSecondOrder(N);
    for (auto const &p: rp.params)
        if (std::isalpha(p.name[0])) model->set_param(p.name, p.value);

    // rp.monitor = new tvb::BoldTVB(N, 720.0, rp.dt, {0});
    // rp.monitor = new tvb::BoldBalloonWindkessel(N, 1.0, 720.0, rp.dt, {0});
    // rp.monitor = new tvb::RawSubSample(1.0, rp.dt, {3});

    // auto *model = new tvb:std:ZerlautAdaptationFirstOrder(N);
    TArray1d sigmas = TArray1d::Constant(model->n_vars(), 0.0);
    for (auto const &p: rp.params)
        if (p.name[0] == '_') {
            auto idx = std::stoi(p.name.substr(2, 1));
            sigmas[idx] = p.value;
        }


    // sigmas << 3e-5, 3e-5, 0.0, 0.0;
    auto *integrator = new tvb::EulerStochastic(rp.dt, new Additive(sigmas, rp.dt));
    // auto *integrator = new tvb::EulerDeterministic();

    auto coupling = new tvb::CouplingLinearSparse(con->weights(), con->delays(), model->cvars());

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

    //auto [step, distance, J_i] = optimize_fic(sim_config, rp.voi, rp.value_base);
    auto [found, a, b, distance, J_i] = optimize_fic_Herzog(sim_config, rp.voi, rp.value_base);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start);

    std::cout << string_format("Optimization time (%s): %d msecs", filename.c_str(), duration.count()) << std::endl;

    total_time += duration;

//    std::vector<long unsigned> shape{J_i.rows(), J_i.cols()};
//    npy::SaveArrayAsNumpy(f_prefix + ".npz", false, shape.size(), shape.data(), (Float*)J_i.data());

    if (found) {

        TArray1d2npz(J_i, f_prefix + ".npz", "J_i");
        TArray1d2npz(TArray1d::Constant(1, 1, model->get_param("G")[0]), f_prefix + ".npz", "G");
        TArray1d2npz(sigmas, f_prefix + ".npz", "s");

        model->set_param("J_i", J_i);
        auto [converged, sim_result] = tvb::simulate(sim_config, 1.0, rp.voi);

        save_fig(sim_result, f_prefix);

        std::ofstream out_txt(f_prefix + ".txt");
        out_txt << string_format("MINIMUM = %f at a=%f, b = %f\n", distance, a, b);


        delete sim_result;
        rp.file_out = filename;
    }

    delete model;
    delete coupling;

    return rp;
}

void run_seq(RunParams rp) {

    tvb::TArray2d C;
    if (rp.file_weights.ends_with(".csv"))
        C = tvb::csv_load(rp.file_weights);
    else if (rp.file_weights.ends_with(".npz"))
        C = tvb::npz2Matrixd(rp.file_weights, "C");
    else
        throw std::runtime_error(string_format("Unknown file extension for: %s", rp.file_weights.c_str()));

    int N = C.rows();

    Float k = 0.15 / (C.rowwise().sum().sum() / N);
    // C *= k;
    // tvb::csv_save("sc_d_norm.csv", C);

    tvb::TArray2d tl;
    if (!rp.file_lengths.empty())
        tl = tvb::csv_load(rp.file_lengths);
    else
        tl = tvb::TArray2d::Zero(C.rows(), C.cols());
    auto *con = new tvb::Connectivity(C, tl, rp.speed);

    TArray1d last_Ji = TArray1d::Constant(N, 1.0);
    float a = 1.0;
    float b = 0.5;
    for (auto G: tvb::range(0.5, 5.0, 45)) {
        string f_prefix = (path(rp.path_out) / path(rp.file_prefix)).lexically_normal().string();
        f_prefix += string_format("_G_%.2f", G);
        string filename = f_prefix + ".png";

        auto *model = new tvb::ReducedWongWangExcInh(N);

        if (!rp.force_output && std::filesystem::exists(filename)) {
            std::cout << string_format("File %s already exists", filename.c_str()) << std::endl;
            TArray2dMap pmap = npz2MatrixdMap(f_prefix + ".npz");
            if (pmap.contains("ab")) {
                TArray1d ab = pmap["ab"];
                a = ab[0];
                b = ab[1];
            }
            if (pmap.contains("J_i"))
                last_Ji = pmap["J_i"].col(0);
            continue;
        }

        std::cout << string_format("Starting computation for: %s", filename.c_str()) << std::endl;

        model->set_param("G", G);
        model->set_param("Ji", last_Ji);

        // rp.monitor = new tvb::BoldTVB(N, 720.0, rp.dt, {0});
        // rp.monitor = new tvb::BoldBalloonWindkessel(N, 1.0, 720.0, rp.dt, {0});
        // rp.monitor = new tvb::RawSubSample(1.0, rp.dt, {3});

        // auto *model = new tvb:std:ZerlautAdaptationFirstOrder(N);
        TArray1d sigmas = TArray1d::Constant(model->n_vars(), 0.0);
        for (auto const &p: rp.params)
            if (p.name[0] == '_') {
                auto idx = std::stoi(p.name.substr(2, 1));
                sigmas[idx] = p.value;
            }

        sigmas[0] = 1e-8;
        sigmas[1] = 1e-8;
        sigmas[2] = 0.0;
        sigmas[3] = 0.0;

        auto *integrator = new tvb::EulerStochastic(rp.dt, new Additive(sigmas, rp.dt));
        // auto *integrator = new tvb::EulerDeterministic();

        auto coupling = new tvb::CouplingLinearSparse(con->weights(), con->delays(), model->cvars());

        SimConfig sim_config;

        sim_config.setModel(model);
        sim_config.setConnectivity(con);
        sim_config.setIntegrator(integrator);
        sim_config.setCoupling(coupling);
        sim_config.setIntegrationInterval(rp.t_start, rp.t_end);
        sim_config.setNumIterations(1);
        sim_config.setDeltaIntegration(0.005);

        //auto [step, distance, J_i] = optimize_fic(sim_config, rp.voi, rp.value_base);
        auto [found, a_best, b_best, distance, J_i] = optimize_fic_Herzog(sim_config, rp.voi, rp.value_base, a, b);


        a = a_best;
        b = b_best;
        last_Ji = J_i;
        TArray1d2npz(J_i, f_prefix + ".npz", "J_i");
        TArray1d2npz(TArray1d::Constant(1, 1, model->get_param("G")[0]), f_prefix + ".npz", "G");
        TArray1d2npz(TArray1d::Constant(1, 1, model->get_param("G")[0]), f_prefix + ".npz", "G");
        TArray1d ab(2);
        ab << a, b;
        TArray1d2npz(ab, f_prefix + ".npz", "ab");

        model->set_param("J_i", J_i);
        auto [converged, sim_result] = tvb::simulate(sim_config, 1.0, rp.voi);

        save_fig(sim_result, f_prefix);

        std::ofstream out_txt(f_prefix + ".txt");
        if (converged)
            out_txt << string_format("MINIMUM = %f at a=%f, b = %f\n", distance, a, b);
        else
            out_txt << string_format("NOT CONVERGED, best minimum = %f at a=%f, b = %f\n", distance, a, b);

        delete sim_result;
        rp.file_out = filename;

        delete model;
        delete coupling;
    }
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
                ("length-matrix", value<std::string>(), "Connection lengths matrix matrix")
                ("speed", value<float>()->default_value(1e6), "Signal speed")
                ("time-start", value<float>()->default_value(0.0), "Start of simulation (ms)")
                ("time-end", value<float>()->default_value(20000.0), "End of simulation (ms)")
                ("dt", value<float>()->default_value(0.1), "Integration step (ms)")
                ("force-output", bool_switch()->default_value(false), "Force overwrite output files")
                ("sigmas", value<std::vector<std::string>>()->multitoken(), "Noise sigmas")
                ("var-of-interest", value<int>()->required(), "Variable of interest in the model to explore")
                ("value-base", value<float>()->default_value(base_value), "Point of equilibrium for excitatory intensity output")
                ("out-path", value<std::string>()->required(), "Output path")
                ("out-file-prefix", value<std::string>()->required(), "Output file prefix");

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
                throw std::runtime_error(string_format("Wrong number of sigma parameters, it should be 4, its %i", params.size() - p_size));
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

        if (true) {
            RunParams pc;
            pc.voi = vm["var-of-interest"].as<int>();
            pc.value_base = vm["value-base"].as<float>();
            pc.file_weights = vm["sc-matrix"].as<std::string>();
            pc.t_start = vm["time-start"].as<float>();
            pc.t_end = vm["time-end"].as<float>();
            pc.dt = vm["dt"].as<float>();
            if (vm.count("length-matrix"))
                pc.file_lengths = vm["length-matrix"].as<std::string>();
            pc.file_prefix = vm["out-file-prefix"].as<std::string>();
            pc.path_out = vm["out-path"].as<std::string>();
            pc.speed = vm["speed"].as<float>();
            pc.force_output = vm["force-output"].as<bool>();
            run_seq(pc);
        } else {
            tvb::ThreadPool<RunParams> tp(6);
            tp.start();
            for (auto &pc: param_combs) {
                pc.voi = vm["var-of-interest"].as<int>();
                pc.value_base = vm["value-base"].as<float>();
                pc.file_weights = vm["sc-matrix"].as<std::string>();
                pc.t_start = vm["time-start"].as<float>();
                pc.t_end = vm["time-end"].as<float>();
                pc.dt = vm["dt"].as<float>();
                if (vm.count("length-matrix"))
                    pc.file_lengths = vm["length-matrix"].as<std::string>();
                pc.file_prefix = vm["out-file-prefix"].as<std::string>();
                pc.path_out = vm["out-path"].as<std::string>();
                pc.speed = vm["speed"].as<float>();
                pc.force_output = vm["force-output"].as<bool>();
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
