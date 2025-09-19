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

#include <tvb-cpp/simulator/integrators/heun_stochastic.h>

using namespace tvb;

State HeunStochastic::scheme(const State &state,
                                 const System &dfun,
                                 const TArray2d &coupling,
                                 const TArray1d &local_coupling,
                                 const TArray2d &stimulus) const {
    // Stochastic Heun method implementation
    // Generate noise for this time step
    TArray2d noise = m_noise->generate(state.rows(), state.cols());
    
    // Step 1: Calculate k1 (drift at the beginning of the interval)
    State k1 = dfun(state, coupling, local_coupling);
    k1 += stimulus;
    
    // Step 2: Calculate intermediate state using Euler step with full noise
    State state_temp = state + this->dt() * k1 + noise;
    
    // Step 3: Calculate k2 (drift at the intermediate state)
    State k2 = dfun(state_temp, coupling, local_coupling);
    k2 += stimulus;
    
    // Step 4: Calculate final state using average of drifts plus noise
    // Note: In stochastic Heun, noise is typically only added once (not averaged)
    State state_next = state + this->dt() * 0.5 * (k1 + k2) + noise;
    
    return state_next;
}