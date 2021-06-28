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
#include <simulator/models/reduced_ww_ext.h>
#include <simulator/integrators/euler_stochastic.h>
#include <simulator/history.h>

int main(int /* argc */, char ** /* argv */ ) {
    std::string baseInPath = "Data_Raw/";

// Configure simulation
    tvb::SimConfig sim_config;

    tvb::Matrixd C = tvb::npz2Matrixd(baseInPath + "Human_66.npz", "C");
    int N = C.rows();

    // Generate random tract lenghts
    tvb::Matrixd tl(N, N);
    tvb::generate(tl, []() { return (5.0 * (double) rand() / (RAND_MAX)); });
    tvb::Connectivity con(C, tl, 1e100);
    auto *ho = new tvb::ReducedWongWangExcInh(N);
// ho->G.fill(we);
    tvb::Vectord sigmas(4);
    sigmas << 3e-5, 3e-5, 0.0, 0.0;
    auto *integrator = new tvb::EulerStochastic(new Additive(sigmas, 0.1));

    sim_config.setModel(ho);
    sim_config.setIntegrator(integrator);
    sim_config.setConnectivity(&con);
    sim_config.setHistory(new tvb::HistoryNoDelays(con.weights(), con.delays(), {3}));
    sim_config.setIntegrationInterval(0.0, 10000.0);
    sim_config.setTimeDelta(0.1);

    auto start = std::chrono::high_resolution_clock::now();
    tvb::StateTrack sresult = simulate(sim_config);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << string_format("Simulation time: %d msecs", duration.count());

}