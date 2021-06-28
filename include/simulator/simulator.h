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

#ifndef TVB_CPP_SIMULATOR_H
#define TVB_CPP_SIMULATOR_H

#include <memory>

#include <simulator/integrator.h>
#include <simulator/coupling.h>
#include <datatypes/connectivity.h>

using namespace std;

namespace tvb {

    class Simulator {
    protected:
        const Coupling *m_coupling;

    public:
        typedef typename std::unique_ptr<Simulator> UPtr;

        Matrixd _loop_compute_node_coupling(int step, History& history) {
            return m_coupling->couple(step, history);
        }

        Vectord _loop_update_stimulus(int step, const Vectord &stimulus) {
            return tvb::Vectord();
        }

        void _loop_update_history(History &history, int step, int n_reg, const State &state) {
            history.update(step, n_reg, state);
        }


        StateTrack *run(const Model* model,
                        const Connectivity* connectivity,
                        Integrator* integrator,
                        const Coupling* coupling,
                        History* history,
                        double start_time, double end_time, double dt,
                        State *initial_state = NULL,
                        int samplingRate = 1);
    };
}

#endif //TVB_CPP_SIMULATOR_H
