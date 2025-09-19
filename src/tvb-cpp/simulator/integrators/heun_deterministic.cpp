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

#include "heun_deterministic.h"

using namespace tvb;

State HeunDeterministic::scheme(const State &state,
                                 const System &dfun,
                                 const TArray2d &coupling,
                                 const TArray1d &local_coupling,
                                 const TArray2d &stimulus) const {
    // Heun's method (improved Euler method)
    // Step 1: Calculate k1 (slope at the beginning of the interval)
    State k1 = dfun(state, coupling, local_coupling);
    // TVB applies stimulus to the first state variable
    k1 += stimulus;
    
    // Step 2: Calculate intermediate state using Euler step
    State state_temp = state + this->dt() * k1;
    
    // Step 3: Calculate k2 (slope at the end of the interval)
    State k2 = dfun(state_temp, coupling, local_coupling);
    k2 += stimulus;
    
    // Step 4: Calculate final state using average of slopes
    State state_next = state + this->dt() * 0.5 * (k1 + k2);
    
//        if self.state_variable_boundaries is not None:
//        self.bound_state(state_next)
//        if self.clamped_state_variable_values is not None:
//        self.clamp_state(state_next)
    return state_next;
}
