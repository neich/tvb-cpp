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

#include <simulator/integrators/euler_deterministic.h>

using namespace tvb;

State EulerDeterministic::scheme(const State &state,
                                 const System &dfun,
                                 const TArray2d &coupling,
                                 const TArray1d &local_coupling,
                                 const TArray1d &stimulus) {
    State d_state = dfun(state, coupling, local_coupling);
    // TVB applies stimulus to the first state variable
    d_state.col(0) += stimulus;
    State t2 = this->dt() * d_state;
    State state_next = state + t2; // 0.0 should be stimulus
//        if self.state_variable_boundaries is not None:
//        self.bound_state(state_next)
//        if self.clamped_state_variable_values is not None:
//        self.clamp_state(state_next)
    return state_next;
}
