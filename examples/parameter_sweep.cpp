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

#include <tools/npz_tools.h>
#include <simulator/simulate.h>
#include <simulator/monitor.h>
#include <simulator/models/reduced_ww_ext.h>
#include <simulator/integrators/euler_stochastic.h>
#include <tools/threadpool.h>
#include "tools/csv_tools.h"
#include "simulator/integrators/euler_deterministic.h"
#include "simulator/models/zerlaut.h"
#include "simulator/simulator.h"

#include <matplotlibcpp.h>
#include <chrono>
#include <filesystem>
#include <thread>

#include <boost/program_options.hpp>
#include <utility>

using namespace boost::program_options;
using namespace tvb;
using namespace std::chrono;
namespace plt = matplotlibcpp;

struct SweepParam {
    string name;
    std::vector<float> values;
};

struct Parameter {
    string name;
    float value;

    Parameter(string name, float value): name(std::move(name)), value(value) {}
};

struct RunParams {
    tvb::RawSubSample* monitor = nullptr;
    std::vector<Parameter> params;
    string file_out;
    string file_prefix;
    string file_weights;
    string file_lengths;
    float speed = 1e6;

    RunParams() = default;
    RunParams(std::vector<Parameter> params, tvb::RawSubSample* monitor): params(std::move(params)), monitor(monitor) {}
};


void save_fig(RunParams &rp);

RunParams run(RunParams rp) {
    string filename = rp.file_prefix;
    for (auto const &p: rp.params) {
        filename += string_format("_%s_%.2f", p.name.c_str(), p.value);
    }
    filename += ".png";

    if (std::filesystem::exists(filename)) {
        std::cout << string_format("File %s already exists", filename.c_str()) << std::endl;
        delete rp.monitor;
        rp.monitor = nullptr;
        return rp;
    }

    tvb::TArray2d C = tvb::csv_load(rp.file_weights);
    int N = C.rows();

    C = C / C.rowwise().sum().maxCoeff() * 2.0;

    tvb::TArray2d tl = tvb::csv_load(rp.file_lengths);
    tvb::Connectivity con(C, tl, rp.speed);

    milliseconds total_time(0);
    std::cout << string_format("Starting computation for: %s", filename.c_str()) << std::endl;


    // auto *ho = new tvb::ReducedWongWangExcInh(N);
    // ho->G.fill(1.0);
    auto *ho = new tvb::ZerlautAdptationSecondOrder(N);
    for (auto const &p: rp.params)
        ho->set_param(p.name, p.value);

        // auto *ho = new tvb:std:ZerlautAdaptationFirstOrder(N);
    // tvb::TArray1d sigmas(4);
    // sigmas << 3e-5, 3e-5, 0.0, 0.0;
    // auto *integrator = new tvb::EulerStochastic(new Additive(sigmas, 0.1));
    auto *integrator = new tvb::EulerDeterministic();

    tvb::SimConfig sim_config;

    sim_config.setModel(ho);
    sim_config.setIntegrator(integrator);
    sim_config.setMonitor(rp.monitor);
    sim_config.setConnectivity(&con);
    auto coupling = new tvb::CouplingLinearSparse(con.weights(), con.delays(), ho->cvars());
    sim_config.setCoupling(coupling);
    // sim_config.setCoupling(new tvb::CouplingLinearDense(con.weights(), con.delays(), ho->cvars()));
    sim_config.setIntegrationInterval(0.0, 10000.0);
    sim_config.setTimeDelta(0.1);

    auto start = std::chrono::high_resolution_clock::now();
    tvb::Simulator simulator{};
    tvb::StateTrack *stateTrack = simulator.run(sim_config.model(),
                                                sim_config.connectivity(),
                                                sim_config.integrator(),
                                                sim_config.monitor(),
                                                sim_config.coupling(),
                                                sim_config.start_time(), sim_config.end_time(), sim_config.dt(),
                                                nullptr,
                                                sim_config.samplingRate());

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start);

    std::cout << string_format("Simulation time (%s): %d msecs", filename.c_str(), duration.count()) << std::endl;

    total_time += duration;

    delete ho;
    delete coupling;
    delete stateTrack;

    rp.file_out = filename;
    return rp;
}


void save_fig(RunParams &rp) {
    if (rp.monitor == nullptr) return;

    size_t t_max = rp.monitor->getRecords().size();
    int N = rp.monitor->getRecords()[0].rows();
    std::vector<std::vector<Float>> y_plot(N, std::vector<Float>(t_max));
    for (unsigned t = 0; t < t_max; ++t)
        for (unsigned n = 0; n < N; ++n)
            y_plot[n][t] = rp.monitor->getRecords()[t](n, 0);

    // tvb::csv_save("./test_simulationRWW_TVB_CPP.csv", y_plot);

    // Plot line from given x and y data. Color is selected automatically.
    std::vector<Float> ls(rp.monitor->getRecords().size());
    std::iota(ls.begin(), ls.end(), 1.0);
    for (unsigned n = 0; n < N; ++n) {
        plt::plot(ls, y_plot[n]);
    }
    // Plot a red dashed line from given x and y data.
    // plt::plot(x, w,"r--");
    // Plot a line whose name will show up as "log(x)" in the legend.

    string title = "Zerlaut 2nd order: ";
    for (auto const &p: rp.params) {
        title += string_format(" %s_%.2f", p.name.c_str(), p.value);
    }

    plt::title(title);
    plt::ylabel("State variable S_e");
    plt::xlabel("Miliseconds");
    // Save the image (file format is determined by the extension)
    plt::save(rp.file_out, 300);

    plt::clf();

    delete rp.monitor;
}

void to_cout(const std::vector<std::string> &v)
{
    std::copy(v.begin(), v.end(), std::ostream_iterator<std::string>{
            std::cout, "\n"});
}

int main(int argc, char **argv) {
    try {
        options_description desc{"Options"};
        desc.add_options()
                ("help,h", "Help screen")
                ("params", value<std::vector<std::string>>()->multitoken(), "Parameters to sweep");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);
        notify(vm);

        std::vector<SweepParam> params;
        if (vm.count("help"))
            std::cout << desc << '\n';
        else if (vm.count("params")) {
            for (auto &s: vm["params"].as<std::vector<std::string>>()) {
                if (std::isalpha(s[0])) {
                    if (!params.empty() && !(params.back().values.size() == 1 ||  params.back().values.size() == 3))
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

        tvb::ThreadPool<RunParams> tp(6);
        tp.start();
        for (auto &pc: param_combs) {
            pc.monitor = new tvb::RawSubSample(10);
//            pc.file_weights = R"(/mnt/d/Dropbox/work/git/research/neuro/tvb/fast_tvb/step2_create_Docker_container/input/gavg_SC_weights.csv)";
//            pc.file_lengths = R"(/mnt/d/Dropbox/work/git/research/neuro/tvb/fast_tvb/step2_create_Docker_container/input/gavg_SC_distances.csv)";
            pc.file_weights = R"(/mnt/d/Dropbox/work/git/research/neuro/tvb/example/hpc_data/jubrain/HCP_N272_SC_JuBrain_294Regions_with_blank_cerebellum/209834_JuBrain_294Regions_10M_ctx_count.csv)";
            pc.file_lengths = R"(/mnt/d/Dropbox/work/git/research/neuro/tvb/example/hpc_data/jubrain/HCP_N272_SC_JuBrain_294Regions_with_blank_cerebellum/209834_JuBrain_294Regions_10M_ctx_length.csv)";
            pc.file_prefix = R"(./simulation_JB_Z)";
            tp.queue_job([pc] { return run(pc); });
        }

        while (!tp.finished()) {
            std::optional<RunParams> op = tp.get_result();
            if (op.has_value())
                save_fig(op.value());
            sleep(1);
        }

        tp.stop();

    }
    catch (const error &ex) {
        std::cerr << ex.what() << '\n';
    }
}
