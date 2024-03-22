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
#include <tvb-cpp/simulator/factory.h>
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
#include <tvb-cpp/tools/cl_params.h>


using namespace tvb;
using namespace std::chrono;
using namespace boost::program_options;

int main(int argc, const char ** argv) {
    std::vector<std::pair<std::string, Float>> params;

    options_description desc{"Options"};
    desc.add_options()
            ("speed", value<tvb::Float>()->default_value(1e6), "Signal speed")
            ("time-start", value<tvb::Float>()->default_value(0.0), "Start of simulation (ms)")
            ("time-end", value<tvb::Float>()->default_value(10000.0), "End of simulation (ms)")
            ("dt", value<tvb::Float>()->default_value(0.1), "Integration step (ms)")
            ("params-file", value<std::string>(), "NPZ file with simulation parameters")
            ("voi", value<int>()->default_value(0), "Index of variable of interest among state variables to explore")
            ("out-file-prefix", value<std::string>()->default_value("out_sim"), "Output file prefix");

    tvb::CLParser cp;
    cp.init(desc, argc, argv);

    auto param_combinations = cp.get_parameter_combinations();
    if (param_combinations.size() > 1)
        throw std::runtime_error("More than 1 parameter combination, you should run a parameter sweep");

    std::string file_prefix = cp.get_option<std::string>("out-file-prefix");

    int N = cp.get_n_rois();
    std::cout << string_format("Connectivity matrix size: %i", N) << std::endl;

    tvb::Connectivity con(*cp.get_sc_matrix(), *cp.get_length_matrix(), cp.get_option<tvb::Float>("speed"));
    auto *model = tvb::Factory::new_model(cp.get_option<std::string>("model"), N);
    auto dt = cp.get_option<tvb::Float>("dt");
    tvb::TArray1d sigmas = cp.get_noise();
    auto *integrator = new tvb::EulerStochastic(dt, new Additive(sigmas, dt));
    int voi = cp.get_option<int>("voi");
    tvb::Float time_end = cp.get_option<tvb::Float>("time-end");

    if (model->n_vars() != sigmas.size())
        throw std::runtime_error(string_format("Provided noise size (%i) does not mathc with model number of state variables (%i)", sigmas.size(), model->n_vars()));

    tvb::Float G = 1.0;
    if (param_combinations.size() == 1)
        cp.init_from_parameters(model, param_combinations[0]);

    // auto *integrator = new tvb::EulerDeterministic(dt);
    auto *coupling = new tvb::CouplingLinearSparse(con.weights(), con.delays(), model->cvars());
    coupling->setScale(G);

    std::cout << string_format("Starting computation for: %s", file_prefix.c_str()) << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    Simulator simulator{};
    TArray2d initial_state = TArray2d::Zero(cp.get_n_rois(), model->n_vars());
    auto *monitor = new TemporalAverage(con.weights().cols(), 1.0, dt, {voi});

    simulator.run(model,
                  &con,
                  integrator,
                  {monitor},
                  coupling,
                  0, time_end,
                  nullptr,
                  &initial_state);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << string_format("Simulation time: %d msecs", duration.count()) << std::endl;

    TArray2d voi_0 = monitor->voi2Array(voi);
    size_t n_records = monitor->getRecords().size();

    tvb::TArray2dMap map_raw;
    map_raw["t_samples"] = tvb::nrange(0.0, 1.0, n_records);
    map_raw["data"] = voi_0;
    tvb::MatrixdMap2npz(file_prefix + "_RAW.npz", map_raw);

    auto *bold_model_tvb = new BoldTVB(2000.0);
    tvb::TArray2dMap map_bold_tvb;
    auto [t_samples_bold_tvb, data_tvb] = bold_model_tvb->compute_bold(voi_0, 1.0);
    map_bold_tvb["t_samples"] = t_samples_bold_tvb;
    map_bold_tvb["data"] = data_tvb;
    tvb::MatrixdMap2npz(file_prefix + "_BOLD.npz", map_bold_tvb);

    return 0;
}