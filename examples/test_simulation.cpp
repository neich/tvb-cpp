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

#include <tvb-root-cpp/tools/npz_tools.h>
#include <tvb-root-cpp/simulator/simulate.h>
#include <tvb-root-cpp/simulator/monitor.h>
#include <tvb-root-cpp/simulator/models/reduced_ww_ext.h>
#include <tvb-root-cpp/simulator/models/montbrio.h>
#include <tvb-root-cpp/simulator/integrators/euler_stochastic.h>
#include <tvb-root-cpp/tools/csv_tools.h>
#include <tvb-root-cpp/simulator/integrators/euler_deterministic.h>
#include <tvb-root-cpp/simulator/models/zerlaut.h>
#include <tvb-root-cpp/simulator/simulator.h>
#include <tvb-root-cpp/simulator/monitors/bold_tvb.h>
#include <tvb-root-cpp/simulator/monitors/bold_BalloonWindkessel.h>

#include <tvb-root-cpp/matplotlibcpp.h>
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
                ("params", value<std::vector<std::string>>()->multitoken(), "Model parameters")
                ("sc-matrix", value<std::string>()->required(), "Structural connectivity matrix")
                ("length-matrix", value<std::string>(), "Connection lengths matrix matrix")
                ("speed", value<float>()->default_value(1e6), "Signal speed")
                ("time-start", value<float>()->default_value(0.0), "Start of simulation (ms)")
                ("time-end", value<float>()->default_value(10000.0), "End of simulation (ms)")
                ("dt", value<float>()->default_value(0.1), "Integration step (ms)")
                ("force-output", bool_switch()->default_value(false), "Force overwrite output files")
                ("params-file", value<std::string>(), "NPZ file with simulation parameters")
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
    string file_weights = vm["sc-matrix"].as<string>();
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
    if (vm.count("length-matrix") > 0)
        tl = tvb::csv_load(vm["length-matrix"].as<string>());
    else
        tl = tvb::TArray2d::Zero(C.rows(), C.cols());
    tvb::Connectivity con(C, tl, vm["speed"].as<float>());

    milliseconds total_time(0);
    std::cout << string_format("Starting computation for: %s", filename.c_str()) << std::endl;

    //auto *model = new tvb::Montbrio(N, rp.t_start, rp.t_end, rp.dt);
    auto *model = new tvb::ReducedWongWangExcInh(N);
    for (auto const &p: params)
        model->set_param_fill(p.first, p.second);

    if (vm.count("params-file") > 0) {
        TArray2dMap pmap = npz2MatrixdMap(vm["params-file"].as<string>());
        if (pmap.contains("G"))
            model->set_param_fill("G", pmap["G"](0, 0));
        if (pmap.contains("J_i"))
            model->set_param_value("J_i", pmap["J_i"].col(0));
    }

    float dt = vm["dt"].as<float>();

    // auto *model = new tvb:std:ZerlautAdaptationFirstOrder(N);
    // tvb::TArray1d sigmas(4);
    // sigmas << 3e-5, 3e-5, 0.0, 0.0;
    // auto *integrator = new tvb::EulerStochastic(new Additive(sigmas, 0.1));
    auto *integrator = new tvb::EulerDeterministic();
    auto *coupling = new tvb::CouplingLinearSparse(con.weights(), con.delays(), model->cvars());

    auto start = std::chrono::high_resolution_clock::now();
    SimConfig sim_config;

    sim_config.setModel(model);
    sim_config.setConnectivity(&con);
    sim_config.setIntegrator(integrator);
    sim_config.setCoupling(coupling);
    sim_config.setIntegrationInterval(vm["time-start"].as<float>(), vm["time-end"].as<float>());
    sim_config.setNumIterations(1);
    sim_config.setDeltaIntegration(0.00001);

    auto [converged, monitor] = tvb::simulate(sim_config, 1.0, 3);

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

    plt::title("RWW");
    plt::ylabel("Ie");
    plt::ylim(0.35, 0.4);
    plt::xlabel("Seconds");
    // Save the image (file format is determined by the extension)
    plt::save("./test_stroke.png", 300);
}