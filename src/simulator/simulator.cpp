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

#include <algorithm>

#include <simulator/simulator.h>

using namespace tvb;

State Simulator::run(Model *model,
                     const Connectivity *connectivity,
                     Integrator *integrator,
                     std::vector<Monitor*> monitors,
                     Coupling *coupling,
                     float start_time, float end_time,
                     State *initial_state) {

    m_coupling = coupling;
    float dt = integrator->dt();

    integrator->configure(start_time, end_time, dt);
    model->init_dependant();

    for (auto mp: monitors)
        mp->setStartTime(start_time);

    auto n_steps = int((end_time - start_time) / dt);

    State state = initial_state != NULL ? *initial_state : model->initial();
    double t = start_time;

    coupling->init(dt, state);

    auto n_reg = connectivity->weights().rows();

    TArray1d local_coupling = TArray1d::Zero(n_reg);
    TArray1d stimulus = TArray1d::Zero(n_reg);

    for (int step = 1; step <= n_steps; ++step) {
        TArray2d node_coupling = this->_loop_compute_node_coupling(step);
        this->_loop_update_stimulus(step, TArray1d()); // TODO: handle stimulus
        state = integrator->scheme(state, *model, node_coupling, local_coupling, stimulus);
#ifndef NDEBUG
        if (tvb::isnan(state))
            throw std::runtime_error("NaN found in integration state!");
#endif
        for (auto mp: monitors)
            mp->record(step, state);
        t += integrator->dt();
        this->_loop_update_history(*coupling, step, state);
    }

    return state;
}


