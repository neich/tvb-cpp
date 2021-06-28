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

#include <simulator/simulator.h>

using namespace tvb;

StateTrack* Simulator::run(const Model *model,
                           const Connectivity *connectivity,
                           Integrator *integrator,
                           const Coupling *coupling,
                           History *history,
                           double start_time, double end_time, double dt,
                           State *initial_state,
                           int samplingRate) {

    m_coupling = coupling;

    integrator->configure(start_time, end_time, dt);

    auto n_steps = int((end_time - start_time) / dt);

    State state = initial_state != NULL ? *initial_state : model->initial();
    StateTrack *result = model->create_track();
    double t = start_time;

    // HistoryDense history(connectivity->weights(), connectivity->delays(), model->cvars());
    history->init(dt, state);

    auto n_reg = connectivity->weights().rows();
    // result->push(state, t);

    Vectord local_coupling = Vectord::Zero(n_reg);
    Vectord stimulus = Vectord::Zero(n_reg);

    for (int step = 1; step <= n_steps; ++step) {
        Matrixd node_coupling = this->_loop_compute_node_coupling(step, *history);
        this->_loop_update_stimulus(step, Vectord()); // TODO: handle stimulus
        state = integrator->scheme(state, *model, node_coupling, local_coupling, stimulus);
        t += integrator->dt();
        if (step % samplingRate == 0)
            result->push(state, t);
        this->_loop_update_history(*history, step, n_reg, state);
        // Vectord output = this->_loop_monitor_output(step, state)

    }

    return result;
}


