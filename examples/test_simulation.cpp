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
#include <tvb-cpp/tools/csv_tools.h>
#include <tvb-cpp/simulator/integrators/euler_deterministic.h>
#include <tvb-cpp/simulator/models/zerlaut.h>
#include <tvb-cpp/simulator/simulator.h>
#include <tvb-cpp/simulator/bold/bold_tvb.h>
#include "tvb-cpp/simulator/bold/bold_BalloonWindkessel.h"
#include "tvb-cpp/simulator/bold/bold_Stephan2007.h"

#include <chrono>
#include <boost/program_options.hpp>
#include <utility>
#include <filesystem>


using namespace tvb;
using namespace std::chrono;
using namespace boost::program_options;

int main(int argc, char ** argv) {


/*
    auto *b = new BoldTVB(2000.0);
    int t_length = 1000000;
    int n_rois = 76;
    tvb::TArray2d tmp = tvb::TArray2d::Zero(t_length, n_rois);
    int val = 0;
    for (int t = 0; t < t_length; ++t)
        for (int n = 0; n < n_rois; ++n) {
            tmp(t, n) = val % 13;
            val++;
        }

    tmp = tvb::npy2Matrixd("test_CNT.npy");
    auto [t_tmp, data_tmp] = b->compute_bold(tmp, 1);

    tvb::Matrixd2np(data_tmp, "tmp_bold_tvb_cpp.npy");
*/

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

    std::string filename = vm["out-file-prefix"].as<std::string>();
    for (auto const &p: params) {
        filename += string_format("_%s_%.2f", p.first.c_str(), p.second);
    }
    filename += ".png";


    if (vm.count("force-output") == 0 && std::filesystem::exists(filename)) {
        std::cout << string_format("File %s already exists", filename.c_str()) << std::endl;
        return 0;
    }

    tvb::TArray2d C;
    std::string file_weights = vm["sc-matrix"].as<std::string>();
    if (file_weights.ends_with(".csv"))
        C = tvb::csv_load(file_weights);
    else if (file_weights.ends_with(".npz"))
        C = tvb::npz2Matrixd(file_weights, "SC");
    else if (file_weights.ends_with(".npy"))
        C = tvb::npy2Matrixd(file_weights);
    else
        throw std::runtime_error(string_format("Unknown file extension for: %s", file_weights.c_str()));

    int N = C.rows();

    C = C / C.maxCoeff() * 0.2;
    // tvb::csv_save("sc_d_norm.csv", C);

    tvb::TArray2d tl;
    if (vm.count("length-matrix") > 0)
        tl = tvb::csv_load(vm["length-matrix"].as<std::string>());
    else
        tl = tvb::TArray2d::Zero(C.rows(), C.cols());
    tvb::Connectivity con(C, tl, vm["speed"].as<float>());

    milliseconds total_time(0);
    std::cout << string_format("Starting computation for: %s", filename.c_str()) << std::endl;

    //auto *model = new tvb::Montbrio(N);
    // tvb::TArray1d sigmas(6);
    // sigmas << 0,0,0,0,1e-3,1e-3;

    // auto *model = new tvb::ReducedWongWangExcInh(N);
    // tvb::TArray1d sigmas(2);
    // sigmas << 1e-5,1e-5;
    auto *model = new tvb::ZerlautAdaptationSecondOrder(N);
    tvb::TArray1d sigmas(8);
    sigmas << 0,0,0,0,0,0,0,1e-5;
    Float G = 1.0;
    auto g_it = std::find_if(params.begin(), params.end(), [](const std::pair<std::string, Float>&p) { return p.first == "G"; });
    if (g_it != params.end()) {
        G = g_it->second;
    }

    for (auto const &p: params)
        if (std::isalpha(p.first[0]) && p.first != "G") model->set_param(p.first, p.second);

    if (vm.count("params-file") > 0) {
        TArray2dMap pmap = npz2MatrixdMap(vm["params-file"].as<std::string>());
        if (pmap.contains("G"))
            G = pmap["G"](0, 0);
        if (pmap.contains("J_i"))
            model->set_param("J_i", pmap["J_i"].col(0));
    }

    model->init_dependant();

    float dt = vm["dt"].as<float>();

    auto *integrator = new tvb::EulerStochastic(dt, new Additive(sigmas, dt));
    // auto *integrator = new tvb::EulerDeterministic(dt);
    auto *coupling = new tvb::CouplingLinearSparse(con.weights(), con.delays(), model->cvars());
    coupling->setScale(G);
    int voi = 0;

    auto start = std::chrono::high_resolution_clock::now();

    Simulator simulator{};
    TArray2d initial_state = TArray2d::Zero(C.cols(), model->n_vars());
    auto *monitor = new TemporalAverage(con.weights().cols(), 1.0 , dt, {0});

    simulator.run(model,
                  &con,
                  integrator,
                  {},
                  coupling,
                  0, 10000,
                  nullptr,
                  &initial_state);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start);

    std::cout << string_format("Simulation time (%s): %d msecs", filename.c_str(), duration.count()) << std::endl;

//    TArray2d voi_0 = monitor->voi2Array(voi);
//    size_t n_records = monitor->getRecords().size();
//
//    auto *bold_model = new BoldStephan2007b(2500.0);
//    auto [t_samples, data] = bold_model->compute_bold(voi_0, 1.0);
//
//    tvb::TArray2dMap map;
//    map["t_samples"] = t_samples;
//    map["data"] = data;
//    tvb::MatrixdMap2npz("test_sim_bold_2007b.npz", map);
//
//    tvb::TArray2dMap map_raw;
//    map_raw["t_samples"] = tvb::nrange(0.0, 1.0, n_records);
//    map_raw["data"] = voi_0;
//    tvb::MatrixdMap2npz("test_sim_raw.npz", map_raw);
//
//    auto *bold_model_tvb = new BoldTVB(2500.0);
//    tvb::TArray2dMap map_bold_tvb;
//    auto [t_samples_bold_tvb, data_tvb] = bold_model_tvb->compute_bold(voi_0, 1.0);
//    map_bold_tvb["t_samples"] = tvb::nrange(0.0, 1.0, n_records);
//    map_bold_tvb["data"] = data_tvb;
//    tvb::MatrixdMap2npz("test_sim_bold_tvb.npz", map_bold_tvb);

    return 0;

    // Plot line from given x and y data. Color is selected automatically.
}