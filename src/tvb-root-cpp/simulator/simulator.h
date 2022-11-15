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

#include <tvb-root-cpp/simulator/integrator.h>
#include <tvb-root-cpp/simulator/coupling.h>
#include <tvb-root-cpp/simulator/stimulus.h>
#include <tvb-root-cpp/simulator/monitor.h>
#include <tvb-root-cpp/datatypes/connectivity.h>

using namespace std;

namespace tvb {

    class Simulator {
    protected:
        const Coupling *m_coupling;

    public:
        typedef typename std::unique_ptr<Simulator> UPtr;

        TArray2d _loop_compute_node_coupling(int step) {
            return m_coupling->couple(step);
        }

        TArray1d _loop_update_stimulus(int step, const TArray1d &stimulus) {
            return tvb::TArray1d();
        }

        void _loop_update_history(Coupling &coupling, int step, const State &state) {
            coupling.update(step, state);
        }


        State run(Model *model,
                  const Connectivity *connectivity,
                  Integrator *integrator,
                  std::vector<Monitor*> monitors,
                  Coupling *coupling,
                  float start_time, float end_time,
                  Stimulus* stimulus = nullptr,
                  State *initial_state = nullptr);
    };
}

#endif //TVB_CPP_SIMULATOR_H
