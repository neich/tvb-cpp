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

StateTrack tvb::simulate(SimConfig &sim_config) {

    // TODO: use monitor to track
//    double max_diff = 1e6;
//
//    Simulator simulator;
//
//    for (unsigned i = 0; i < sim_config.num_iterations(); ++i) {
//        simulator.run(sim_config.model(),
//                       sim_config.connectivity(),
//                       sim_config.integrator(),
//                       sim_config.monitor(),
//                       sim_config.coupling(),
//                       sim_config.start_time(), sim_config.end_time(), sim_config.dt(),
//                       NULL,
//                       sim_config.samplingRate());
//
//        double v_max = -1e10;
//        double v_min = 1e10;
//        for (unsigned i = stateTrack->states().size() - 1000; i < stateTrack->states().size(); ++i) {
//            double s_max = stateTrack->states()[i].col(sim_config.svar_index()).maxCoeff();
//            double s_min = stateTrack->states()[i].col(sim_config.svar_index()).minCoeff();
//            v_max = std::max(v_max, s_max);
//            v_min = std::min(v_min, s_min);
//        }
//
//        double v_diff = v_max - v_min;
//        max_diff = std::max(max_diff, v_diff);
//
//        full_track.append(*stateTrack);
//
//        if (max_diff < sim_config.delta_integration())
//            break;
//
//    }
//
//    StateTrack result;
//    result.m_states = full_track.states();
//    result.m_times = full_track.times();
//    // result.m_max_diff = max_diff;
//
//    return result;
    return StateTrack();
}
