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
#include <tvb-cpp/simulator/monitors/bold_tvb.h>
#include <tvb-cpp/simulator/monitors/bold_BalloonWindkessel.h>
#include "zerlaut_gaba.h"

#include <chrono>
#include <boost/program_options.hpp>
#include <utility>
#include <filesystem>


using namespace tvb;
using namespace std::chrono;
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
                ("gaba-vector", value<std::string>(), "Vector with neuroreceptor density")
               // ("norm", value<std::vector<std::string>>()->multitoken(), "Matrix normalilzation method")
                ("out-file-prefix", value<std::string>()->required(), "Output file prefix");

        store(command_line_parser(argc, argv)
                      .options(desc)
                      .style(command_line_style::unix_style ^ command_line_style::allow_short)
                      .run(), vm);
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

    tvb::Float G = 1.0;
    auto g_it = std::find_if(params.begin(), params.end(), [](const std::pair<std::string, Float>& p) { return p.first == "G"; });
    if (g_it != params.end()) {
        G = g_it->second;
        params.erase(g_it);
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

    double maxC = C.rowwise().sum().maxCoeff();
    C = C / maxC;

    tvb::TArray2d tl;
    if (vm.count("length-matrix") > 0)
        tl = tvb::csv_load(vm["length-matrix"].as<std::string>());
    else
        tl = tvb::TArray2d::Zero(C.rows(), C.cols());
    tvb::Connectivity con(C, tl, vm["speed"].as<float>());


    //auto *model = new tvb::Montbrio(N, rp.t_start, rp.t_end, rp.dt);
    // auto *model = new tvb::ReducedWongWangExcInh(N);
    auto *model = new ZerlautGABA(N);
    model->configure();
    for (auto const &p: params)
        model->set_param(p.first, p.second);

    TArray1d gaba_vector;
    std::string file_gaba = vm["gaba-vector"].as<std::string>();
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

    model->set_param("gaba_ratio", gaba_vector);
    model->init_dependant();

//    for (auto &s: model->get_param_list())
//        std::cout << string_format("Parameter %s = %f\n", s.c_str(), model->get_param_value(s)[0]);

    // std::cout << std::flush;

    if (vm.count("params-file") > 0) {
        TArray2dMap pmap = npz2MatrixdMap(vm["params-file"].as<std::string>());
        if (pmap.contains("G"))
            model->set_param("G", pmap["G"](0, 0));
        if (pmap.contains("J_i"))
            model->set_param("J_i", pmap["J_i"].col(0));
    }

    float dt = vm["dt"].as<float>();

    // auto *model = new tvb:std:ZerlautAdaptationFirstOrder(N);
    tvb::TArray1d sigmas(8);
    sigmas << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1;
    auto *integrator = new tvb::EulerStochastic(dt, new Additive(sigmas, dt));
    // auto *integrator = new tvb::EulerDeterministic(dt);
    auto *coupling = new tvb::CouplingLinearSparse(con.weights(), con.delays(), model->cvars());
    coupling->setScale(G);

    auto start = std::chrono::high_resolution_clock::now();
    SimConfig sim_config;

    sim_config.setModel(model);
    sim_config.setConnectivity(&con);
    sim_config.setIntegrator(integrator);
    sim_config.setCoupling(coupling);
    sim_config.setIntegrationInterval(vm["time-start"].as<float>(), vm["time-end"].as<float>());
    sim_config.setNumIterations(1);
    sim_config.setDeltaIntegration(0.00001);

    Simulator simulator{};
    BoldTVB *btvb = new BoldTVB(con.weights().cols(), 2500, dt, {0});
    TemporalAverage* ta_monitor = new TemporalAverage(con.weights().cols(), 250.0, dt, {0});
    TArray2d initial_state = TArray2d::Zero(C.cols(), model->n_vars());

    simulator.run(sim_config.model(),
                  sim_config.connectivity(),
                  sim_config.integrator(),
                  {btvb, ta_monitor},
                  sim_config.coupling(),
                  0, vm["time-end"].as<float>(),
                  nullptr,
                  &initial_state);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start);

    std::cout << string_format("Simulation time (%s): %d msecs", filename.c_str(), duration.count()) << std::endl;

    TArray2d bold_data = btvb->voi2Array(0);
    TArray2d raw_data = ta_monitor->voi2Array(0);
    std::cout << string_format("Bold data with %d ROIs an %d samples", bold_data.rows(), bold_data.cols());

    Matrixd2np(bold_data.transpose(), "bold_data.npy");
    Matrixd2np(raw_data.transpose(), "raw_data.npy");

}