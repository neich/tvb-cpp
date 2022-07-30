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

using namespace tvb;
using namespace std::chrono;
namespace plt = matplotlibcpp;

typedef std::tuple<float, float, tvb::RawSubSample*> Params;

void save_fig(Params &tuple);

Params run(Params params) {
    auto [qe, qi, monitor] = params;
    string filename = string_format("./test_simulationZ_qe_%.2f_qi_%.2f.png", qe, qi);
    if (std::filesystem::exists(filename)) {
        std::cout << string_format("Simulation (qe: %.2f, qi: %.2f) already done", qe, qi) << std::endl;
        delete monitor;
        return {qe, qi, NULL};
    }

    tvb::TArray2d C = tvb::csv_load(R"(/mnt/d/Dropbox/work/git/research/neuro/tvb/fast_tvb/step2_create_Docker_container/input/gavg_SC_weights.csv)");
    int N = C.rows();

    // C = C / C.rowwise().sum().maxCoeff() * 2.0;

    tvb::TArray2d tl = tvb::csv_load(R"(/mnt/d/Dropbox/work/git/research/neuro/tvb/fast_tvb/step2_create_Docker_container/input/gavg_SC_distances.csv)");

    tvb::Connectivity con(C, tl, 12.5);

    milliseconds total_time(0);


    std::cout << string_format("Starting (qe: %.2f, qi: %.2f)", qe, qi) << std::endl;


    // auto *ho = new tvb::ReducedWongWangExcInh(N);
    // ho->G.fill(1.0);
    auto *ho = new tvb::ZerlautAdptationSecondOrder(N);
    ho->Q_i.fill(qi);
    ho->Q_e.fill(qe);
    // auto *ho = new tvb:std:ZerlautAdaptationFirstOrder(N);
    // tvb::TArray1d sigmas(4);
    // sigmas << 3e-5, 3e-5, 0.0, 0.0;
    // auto *integrator = new tvb::EulerStochastic(new Additive(sigmas, 0.1));
    auto *integrator = new tvb::EulerDeterministic();

    tvb::SimConfig sim_config;

    sim_config.setModel(ho);
    sim_config.setIntegrator(integrator);
    sim_config.setMonitor(monitor);
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
                                                NULL,
                                                sim_config.samplingRate());

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start);

    std::cout << string_format("Simulation time (qe: %.2f, qi: %.2f): %d msecs", duration.count(), qe, qi) << std::endl;

    total_time += duration;

//    size_t t_max = monitor->getRecords().size();
//    std::vector<std::vector<Float>> y_plot(N, std::vector<Float>(t_max));
//    for (unsigned t = 0; t < t_max; ++t)
//        for (unsigned n = 0; n < N; ++n)
//            y_plot[n][t] = monitor->getRecords()[t](n, 0);
//
//    // tvb::csv_save("./test_simulationRWW_TVB_CPP.csv", y_plot);
//
//    // Plot line from given x and y data. Color is selected automatically.
//    std::vector<Float> ls(monitor->getRecords().size());
//    std::iota(ls.begin(), ls.end(), 1.0);
//    for (unsigned n = 0; n < N; ++n) {
//        plt::plot(ls, y_plot[n]);
//    }
//    // Plot a red dashed line from given x and y data.
//    // plt::plot(x, w,"r--");
//    // Plot a line whose name will show up as "log(x)" in the legend.
//
//    plt::title(string_format("Zerlaut 2nd order: Qe: %.2f, Qi: %.2f", qe, qi));
//    plt::ylabel("State variable S_e");
//    plt::xlabel("Miliseconds");
//    // Save the image (file format is determined by the extension)
//    plt::save(filename, 300);

    delete ho;
    delete coupling;
    delete stateTrack;

    return params;
}


void save_fig(Params &params) {
    auto [qe, qi, monitor] = params;

    if (monitor == NULL) return;

    string filename = string_format("./test_simulationZ_qe_%.2f_qi_%.2f.png", qe, qi);

    size_t t_max = monitor->getRecords().size();
    int N = monitor->getRecords()[0].rows();
    std::vector<std::vector<Float>> y_plot(N, std::vector<Float>(t_max));
    for (unsigned t = 0; t < t_max; ++t)
        for (unsigned n = 0; n < N; ++n)
            y_plot[n][t] = monitor->getRecords()[t](n, 0);

    // tvb::csv_save("./test_simulationRWW_TVB_CPP.csv", y_plot);

    // Plot line from given x and y data. Color is selected automatically.
    std::vector<Float> ls(monitor->getRecords().size());
    std::iota(ls.begin(), ls.end(), 1.0);
    for (unsigned n = 0; n < N; ++n) {
        plt::plot(ls, y_plot[n]);
    }
    // Plot a red dashed line from given x and y data.
    // plt::plot(x, w,"r--");
    // Plot a line whose name will show up as "log(x)" in the legend.

    plt::title(string_format("Zerlaut 2nd order: Qe: %.2f, Qi: %.2f", qe, qi));
    plt::ylabel("State variable S_e");
    plt::xlabel("Miliseconds");
    // Save the image (file format is determined by the extension)
    plt::save(filename, 300);

    plt::clf();
}

int main(int /* argc */, char ** /* argv */ ) {

    tvb::ThreadPool<Params> tp(6);
    tp.start();
    std::vector<tvb::RawSubSample*> monitors;
    for (auto qe: tvb::range(1.0f, 3.0f, 8)) {
        for (auto qi: tvb::range(2.0f, 6.0f, 8)) {
            std::cout << string_format("Queueing job (qe: %.2f, qi: %.2f)", qe, qi) << std::endl;
            auto *monitor = new tvb::RawSubSample(10);
            tp.queue_job([qe, qi, monitor] { return run({qe, qi, monitor}); });
        }
    }
    while (!tp.empty()) {
        std::optional<Params> op = tp.get_result();
        if (op.has_value())
            save_fig(op.value());
        sleep(1);
    }
    while (tp.has_results()) {
        std::optional<Params> op = tp.get_result();
        if (op.has_value())
            save_fig(op.value());
    }

    tp.stop();
}
