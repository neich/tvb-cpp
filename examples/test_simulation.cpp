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
#include <simulator/models/montbrio.h>
#include <simulator/integrators/euler_stochastic.h>
#include "tools/csv_tools.h"
#include "simulator/integrators/euler_deterministic.h"
#include "simulator/models/zerlaut.h"
#include "simulator/simulator.h"
#include "simulator/monitors/bold_tvb.h"
#include "simulator/monitors/bold_BalloonWindkessel.h"

#include <matplotlibcpp.h>
#include <chrono>
#include <boost/program_options.hpp>
#include <utility>
#include <filesystem>


using namespace tvb;
using namespace std::chrono;
namespace plt = matplotlibcpp;
using namespace boost::program_options;

int main(int argc, char ** argv) {

    std::vector<std::pair<std::string, Float>> params;
    variables_map vm;
    try {
        options_description desc{"Options"};
        desc.add_options()
                ("help,h", "Help screen")
                ("params", value<std::vector<std::string>>()->multitoken()->required(), "Parameters to sweep")
                ("sc-matrix", value<std::string>()->required(), "Structural connectivity matrix")
                ("length-matrix", value<std::string>(), "Connection lengths matrix matrix")
                ("speed", value<float>()->default_value(1e6), "Signal speed")
                ("time-start", value<float>()->default_value(0.0), "Start of simulation (ms)")
                ("time-end", value<float>()->default_value(10000.0), "End of simulation (ms)")
                ("dt", value<float>()->default_value(0.1), "Integration step (ms)")
                ("force-output", bool_switch()->default_value(false), "Force overwrite output files")
                ("out-file-prefix", value<std::string>()->required(), "Output file prefix");

        store(parse_command_line(argc, argv, desc), vm);
        notify(vm);

        if (vm.count("help"))
            std::cout << desc << '\n';
        else if (vm.count("params")) {
            bool param_started = false;
            for (auto &s: vm["params"].as<std::vector<std::string>>()) {
                if (std::isalpha(s[0])) {
                    if (param_started)
                        throw std::runtime_error(string_format("Invalid value for parameter %s", params.back().first.c_str()));
                    params.emplace_back();
                    params.back().first = s;
                    param_started = true;
                } else {
                    if (!param_started)
                        throw std::runtime_error(string_format("Value %f has no parameter associated", s.c_str()));
                    try {
                        float value = std::stof(s);
                        params.back().second = value;
                        param_started = false;
                    } catch (const error &ex) {
                        throw std::runtime_error("Syntax error in sweep parameters");
                    }
                }
            }
            if (param_started)
                throw std::runtime_error(string_format("Parameter %s has no associated value\n", params.back().first.c_str()));
        }
    }
    catch (const error &ex) {
        std::cerr << ex.what() << '\n';
    }

    string filename = vm["out-file-prefix"].as<string>();
    for (auto const &p: params) {
        filename += string_format("_%s_%.2f", p.first.c_str(), p.second);
    }
    filename += ".png";


    if (vm.count("force-output") == 0 && std::filesystem::exists(filename)) {
        std::cout << string_format("File %s already exists", filename.c_str()) << std::endl;
        return 0;
    }

    tvb::TArray2d C;
    string file_weights = vm["file-weights"].as<string>();
    if (file_weights.ends_with(".csv"))
        C = tvb::csv_load(file_weights);
    else if (file_weights.ends_with(".npz"))
        C = tvb::npz2Matrixd(file_weights, "SC");
    else
        throw std::runtime_error(string_format("Unknown file extension for: %s", file_weights.c_str()));

    int N = C.rows();

    // C = C / C.rowwise().sum().maxCoeff() * 2.0;
    // tvb::csv_save("sc_d_norm.csv", C);

    tvb::TArray2d tl;
    string file_lengths = vm["file_lengths"].as<string>();
    if (!file_lengths.empty())
        tl = tvb::csv_load(file_lengths);
    else
        tl = tvb::TArray2d::Zero(C.rows(), C.cols());
    tvb::Connectivity con(C, tl, vm["speed"].as<float>());

    milliseconds total_time(0);
    std::cout << string_format("Starting computation for: %s", filename.c_str()) << std::endl;

    //auto *model = new tvb::Montbrio(N, rp.t_start, rp.t_end, rp.dt);
    auto *model = new tvb::ReducedWongWangExcInh(N);
    // auto *model = new tvb::ZerlautAdptationSecondOrder(N);
    for (auto const &p: params)
        model->set_param(p.first, p.second);

    float dt = vm["dt"].as<float>();

    auto *monitor = new tvb::BoldTVB(N, 720.0, dt, {0});
    // rp.monitor = new tvb::BoldBalloonWindkessel(N, 1.0, 720.0, rp.dt, {0});
    // rp.monitor = new tvb::RawSubSample(1.0, rp.dt, {0});

    // auto *model = new tvb:std:ZerlautAdaptationFirstOrder(N);
    // tvb::TArray1d sigmas(4);
    // sigmas << 3e-5, 3e-5, 0.0, 0.0;
    // auto *integrator = new tvb::EulerStochastic(new Additive(sigmas, 0.1));
    auto *integrator = new tvb::EulerDeterministic();
    auto *coupling = new tvb::CouplingLinearSparse(con.weights(), con.delays(), model->cvars());

    auto start = std::chrono::high_resolution_clock::now();
    tvb::Simulator simulator{};
    simulator.run(model,
                  &con,
                  integrator,
                  monitor,
                  coupling,
                  vm["time-start"].as<float>(), vm["time-end"].as<float>(), dt,
                  nullptr);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start);

    std::cout << string_format("Simulation time (%s): %d msecs", filename.c_str(), duration.count()) << std::endl;


    size_t n_records = monitor->getRecords().size();
    std::vector<std::vector<Float>> y_plot(N, std::vector<Float>(n_records));
    for (unsigned t = 0; t < n_records; ++t)
        for (unsigned n = 0; n < N; ++n)
            y_plot[n][t] = monitor->getRecords()[t].record(n, 0);

    // tvb::csv_save("./paper_RWW_BOLD_TVBCPP.csv", y_plot);

    // Plot line from given x and y data. Color is selected automatically.
    std::vector<Float> ls(n_records);
    std::transform(monitor->getRecords().begin(), monitor->getRecords().end(), ls.begin(),
                   [](const Monitor::Record &r) { return r.time/1000; });
    for (unsigned n = 0; n < N; ++n) {
        plt::plot(ls, y_plot[n]);
    }
    // Plot a red dashed line from given x and y data.
    // plt::plot(x, w,"r--");
    // Plot a line whose name will show up as "log(x)" in the legend.

    plt::title("Montbrio");
    plt::ylabel("r_e");
    plt::xlabel("Seconds");
    // Save the image (file format is determined by the extension)
    plt::save("./test_M.png", 300);
}