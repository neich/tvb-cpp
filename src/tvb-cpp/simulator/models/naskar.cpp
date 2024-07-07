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

#include <tvb-cpp/simulator/models/naskar.h>

using namespace tvb;

State Naskar::operator()(const State &x,
        const TArray2d &coupling,
        const TArray1d &local_coupling) const {

    State derivative(m_n_nodes, m_n_vars);

    const TArray1d &S_e = x.col(0);
    const TArray1d &S_i = x.col(1);
    const TArray1d &J = x.col(2);

    TArray1d lc_0 = local_coupling * S_e;

    TArray1d total_coupling = coupling.col(0) + lc_0;

    TArray1d Ie = W_e * I0 + w * J_N * S_e + J_N * total_coupling - J * S_i + I_external;
    TArray1d Ii = W_i * I0 + J_N * S_e - S_i;
    TArray1d y = M_e * (a_e * Ie - b_e);
    TArray1d re = y / (1. - (-d_e * y).exp());
    y = M_i * (a_i * Ii - b_i);
    TArray1d ri = y / (1. - (-d_i * y).exp());

    derivative.col(0) = -S_e * B_e + alpha_e * t_glu * (1. - S_e) * re / 1000.;
    derivative.col(1) = -S_i * B_i + alpha_i * t_gaba * (1. - S_i) * ri / 1000.;
    derivative.col(2) = gamma * ri / 1000. * (re - rho) / 1000.;

#ifndef NDEBUG
    if (!derivative.allFinite())
        throw std::runtime_error("Non finite value found while computing Naskar derivative");
#endif

    return derivative;
}

