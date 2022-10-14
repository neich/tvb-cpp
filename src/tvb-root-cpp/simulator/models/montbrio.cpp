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

#include <tvb-root-cpp/simulator/models/montbrio.h>

using namespace tvb;

State Montbrio::operator()(const State &x,
                           const TArray2d &coupling,
                           const TArray1d &local_coupling) {

    State derivative(m_n_nodes, m_n_vars);


    const TArray1d &r_e = x.col(0);
    const TArray1d &r_i = x.col(1);
    const TArray1d &u_e = x.col(2);
    const TArray1d &u_i = x.col(3);
    const TArray1d &S_ee = x.col(4);
    const TArray1d &S_ei = x.col(5);
    const TArray1d &S_ie = x.col(6);
    const TArray1d &S_ii = x.col(7);
    const TArray1d &c_0 = coupling.col(0);

    derivative.col(0) = (delta_e / (M_PI*tau_e) + Float(2.0) * r_e * u_e) / tau_e; // no coupling
    derivative.col(1) = (delta_i / (M_PI*tau_i) + Float(2.0) * r_i * u_i) / tau_i; // no coupling
    derivative.col(2) = (eta_e + u_e.pow(2.0) - (r_e * M_PI * tau_e).pow(2.0) + I_e) / tau_e;
    derivative.col(3) = (eta_i + u_i.pow(2.0) - (r_i * M_PI * tau_i).pow(2.0) + I_i) / tau_i;
    derivative.col(4) = -S_ee + J_ee*r_e;
    derivative.col(5) = -S_ei + J_ei*r_i;
    derivative.col(6) = -S_ie + J_ie*r_e;
    derivative.col(7) = -S_ii + J_ii*r_i;
    I_e = I_ext + tau_e*S_ee - J_e*S_ei + J_ee * G * c_0;
    I_i = I_ext + tau_i*S_ie - tau_i*S_ii;

    return derivative;

}

