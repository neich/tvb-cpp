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

#include <chrono>

#include <simulator/simulate.h>
#include <simulator/simulator.h>

using namespace tvb;

std::tuple<bool, Monitor*> tvb::simulate(SimConfig &sim_config, float sub_period, int voi) {

    // TODO: use monitor to track
    float max_diff = -1e6;

    Simulator simulator;
    int N = sim_config.connectivity()->weights().rows();

    TArray2d last_state;
    Monitor *monitor = new AverageSubSample(N, sub_period, sim_config.integrator()->dt(), {voi});
    sim_config.addMonitor(monitor);
    bool converged = false;
    float t_start = sim_config.start_time();
    float t_end = sim_config.end_time();
    float t_interval = t_end - t_start;
    for (unsigned i = 0; i < sim_config.num_iterations(); ++i) {
        last_state = simulator.run(sim_config.model(),
                       sim_config.connectivity(),
                       sim_config.integrator(),
                       sim_config.monitors(),
                       sim_config.coupling(),
                       t_start, t_end,
                       i == 0 ? nullptr : &last_state);

        TArray1d v_max = TArray1d::Constant(N, -1e6);
        TArray1d v_min = TArray1d::Constant(N, 1e6);
        auto const& records = monitor->getRecords();
        auto n_records = records.size();
        for (unsigned i = n_records - 4000; i < n_records; ++i) {
            const TArray1d &r = records[i].record.col(0);
            for (int n = 0; n < N; n++) {
                v_max(n) = std::max(v_max(n), r(n));
                v_min(n) = std::min(v_min(n), r(n));
            }
        }


        TArray1d v_diff = v_max - v_min;
        max_diff = std::max(max_diff, v_diff.maxCoeff());

        if (max_diff < sim_config.delta_integration()) {
            converged = true;
            break;
        }

        t_start += t_interval;
        t_end += t_interval;
    }

    sim_config.removeMonitor(monitor);
    return {converged, monitor};
}
