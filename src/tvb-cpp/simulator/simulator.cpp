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

#include <cassert>
#include <algorithm>

#include "simulator.h"

using namespace tvb;

State Simulator::run(const Model *model,
                     const Connectivity *connectivity,
                     const Integrator *integrator,
                     const std::vector<Monitor*>& monitors,
                     Coupling *coupling,
                     float start_time, float end_time,
                     Stimulus* stimulus,
                     State *initial_state) {

    assert(("Model cannot be null!", model != nullptr));
    assert(("Connectivity cannot be null!", connectivity != nullptr));
    assert(("Integrator cannot be null!", integrator != nullptr));
    assert(("There has to be at least a monitor!", monitors.size() > 0));
    assert(("Coupling cannot be null!", coupling != nullptr));

    m_coupling = coupling;
    float dt = integrator->dt();

    if (stimulus) stimulus->configure(start_time, end_time, dt);

    for (auto mp: monitors)
        mp->setStartTime(start_time);

    auto n_steps = int((end_time - start_time) / dt);

    State state = initial_state != nullptr ? *initial_state : model->initial();
    double t = start_time;

    coupling->init(dt, state);

    auto n_reg = connectivity->weights().rows();

    TArray1d local_coupling = TArray1d::Zero(n_reg);
    TArray2d current_st = stimulus != nullptr ? stimulus->initial(state) : TArray2d::Zero(state.rows(), state.cols());

    for (int step = 1; step <= n_steps; ++step) {
        TArray2d node_coupling = this->_loop_compute_node_coupling(step);
        if (stimulus)
            current_st = stimulus->update(step, current_st);
        state = integrator->scheme(state, *model, node_coupling, local_coupling, current_st);
#ifndef NDEBUG
        if (!state.allFinite())
            throw std::runtime_error(string_format("NaN found in integration state, step = %d", step));
#endif
        for (auto mp: monitors)
            mp->record(step, state);
        t += integrator->dt();
        this->_loop_update_history(*coupling, step, state);
    }

    return state;
}