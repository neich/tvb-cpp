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
#include "tools/csv_tools.h"
#include "simulator/integrators/euler_deterministic.h"
#include "simulator/models/zerlaut.h"
#include "simulator/simulator.h"

#include <matplotlibcpp.h>
#include <chrono>

using namespace tvb;
using namespace std::chrono;
namespace plt = matplotlibcpp;

int main(int /* argc */, char ** /* argv */ ) {

    std::string baseInPath = "Data_Raw/";

// Configure simulation
    tvb::SimConfig sim_config;

    // tvb::TArray2d C = tvb::csv_load(R"(D:\Dropbox\work\git\research\neuro\tvb\example\hpc_data\jubrain\HCP_N272_SC_JuBrain_294Regions_with_blank_cerebellum\101309_JuBrain_294Regions_10M_ctx_count.csv)");
    // tvb::TArray2d C = tvb::csv_load(R"(/mnt/d/Dropbox/work/git/research/neuro/tvb/example/hpc_data/jubrain/HCP_N272_SC_JuBrain_294Regions_with_blank_cerebellum/101309_JuBrain_294Regions_10M_ctx_count.csv)");
    tvb::TArray2d C = tvb::csv_load(R"(/mnt/d/Dropbox/work/git/research/neuro/tvb/fast_tvb/step2_create_Docker_container/input/gavg_SC_weights.csv)");
    int N = C.rows();

    // C = C / C.rowwise().sum().maxCoeff() * 2.0;

    tvb::TArray2d tl = tvb::csv_load(R"(/mnt/d/Dropbox/work/git/research/neuro/tvb/fast_tvb/step2_create_Docker_container/input/gavg_SC_distances.csv)");

    milliseconds total_time(0);

    for (int lap = 0; lap < 1; ++lap) {
        // Generate random tract lenghts
        tvb::Connectivity con(C, tl, 10000);
        // auto *ho = new tvb::ReducedWongWangExcInh(N);
        // ho->G.fill(1.0);
        auto *ho = new tvb::ZerlautAdptationSecondOrder(N);
        // auto *ho = new tvb::ZerlautAdaptationFirstOrder(N);
        tvb::TArray1d sigmas(4);
        sigmas << 3e-5, 3e-5, 0.0, 0.0;
        // auto *integrator = new tvb::EulerStochastic(new Additive(sigmas, 0.1));
        auto *integrator = new tvb::EulerDeterministic();

        sim_config.setModel(ho);
        sim_config.setIntegrator(integrator);
        auto *monitor = new tvb::RawSubSample(10);
        sim_config.setMonitor(monitor);
        sim_config.setConnectivity(&con);
        sim_config.setCoupling(new tvb::CouplingLinearSparse(con.weights(), con.delays(), ho->cvars()));
        // sim_config.setCoupling(new tvb::CouplingLinearDense(con.weights(), con.delays(), ho->cvars()));
        sim_config.setIntegrationInterval(0.0, 3000.0);
        sim_config.setTimeDelta(0.1);

        auto start = std::chrono::high_resolution_clock::now();
        tvb::Simulator simulator;
        tvb::StateTrack *stateTrack = simulator.run(sim_config.model(),
                                               sim_config.connectivity(),
                                               sim_config.integrator(),
                                               sim_config.monitor(),
                                               sim_config.coupling(),
                                               sim_config.start_time(), sim_config.end_time(), sim_config.dt(),
                                               NULL,
                                               sim_config.samplingRate());

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - start);

        std::cout << string_format("Simulation time: %d msecs", duration.count()) << std::endl;

        total_time += duration;

        size_t t_max = monitor->getRecords().size();
        std::vector<std::vector<Float>> y_plot(N, std::vector<Float>(t_max));
        for (unsigned t = 0; t < t_max; ++t)
            for (unsigned n = 0; n < N; ++n)
                y_plot[n][t] = monitor->getRecords()[t](n, 0);

        tvb::csv_save("./test_simulationRWW_TVB_CPP.csv", y_plot);

        // Plot line from given x and y data. Color is selected automatically.
        std::vector<Float> ls(monitor->getRecords().size());
        std::iota(ls.begin(), ls.end(), 1.0);
        for (unsigned n = 0; n < N; ++n) {
            plt::plot(ls, y_plot[n]);
        }
        // Plot a red dashed line from given x and y data.
        // plt::plot(x, w,"r--");
        // Plot a line whose name will show up as "log(x)" in the legend.

        plt::title("Reduced Wong Wang - TVB C++");
        plt::ylabel("State variable S_e");
        plt::xlabel("Miliseconds");
        // Save the image (file format is determined by the extension)
        plt::save("./test_simulationZ_TVB_CPP.png", 300);

    }

    std::cout << string_format("Average simulation time: %d msecs", total_time.count()) << std::endl;

//    size_t t_max = monitor->getRecords().size();
//    std::vector<std::vector<double>> y_plot(N, std::vector<double>(t_max));
//    for (unsigned t = 0; t < t_max; ++t)
//        for (unsigned n = 0; n < N; ++n)
//            y_plot[n][t] = monitor->getRecords()[t](n, 0);
//
//    // Plot line from given x and y data. Color is selected automatically.
//    std::vector<double> ls(monitor->getRecords().size());
//    std::iota(ls.begin(), ls.end(), 1.0);
//    for (unsigned n = 0; n < N; ++n) {
//        plt::plot(ls, y_plot[n]);
//    }
//    // Plot a red dashed line from given x and y data.
//    // plt::plot(x, w,"r--");
//    // Plot a line whose name will show up as "log(x)" in the legend.
//
//    plt::title("Reduced Wong Wang - TVB C++");
//    plt::ylabel("State variable S_e");
//    plt::xlabel("Miliseconds");
//    // Save the image (file format is determined by the extension)
//    plt::save("./test_simulationRWW_TVB_CPP.png", 300);

}