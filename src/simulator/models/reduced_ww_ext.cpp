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

#include <simulator/models/reduced_ww_ext.h>

using namespace tvb;

State ReducedWongWangExcInh::operator()(const State &x,
                                        const Matrixd &coupling,
                                        const Vectord &local_coupling) const {

    State derivative(m_n_nodes, m_n_vars);


    const Vectord &S_e = x.col(0);
    const Vectord &S_i = x.col(1);

    Vectord lc_0 = local_coupling * S_e;

    Vectord total_coupling = this->G * this->J_N * (coupling.col(0) + lc_0);

    Vectord J_N_S_e = this->J_N * S_e;

// double I_ext = 0.0;  // Resting state

    Vectord inh = this->J_i * S_i;

    Vectord I_e = this->W_e * this->I_o + this->w_p * J_N_S_e + total_coupling - inh; //  + I_ext;

    Vectord x_e = this->a_e * I_e - this->b_e;
    Vectord tmp_x_e_d = 1.0 - (-this->d_e * x_e).exp();
    Vectord H_e = x_e / tmp_x_e_d;

    derivative.col(0) = -(S_e / this->tau_e) + (1.0 - S_e) * H_e * this->gamma_e;

    Vectord I_i = this->W_i * this->I_o + J_N_S_e - S_i + this->lambda * total_coupling;

    Vectord x_i = this->a_i * I_i - this->b_i;
    Vectord tmp_x_i_d = 1.0 - (-this->d_i * x_i).exp();
    Vectord H_i = x_i / tmp_x_i_d;

    derivative.col(1) = -(S_i / this->tau_i) + H_i * this->gamma_i;

    derivative.col(2) = H_e - x.col(2);
    derivative.col(3) = I_e - x.col(3);

    return derivative;

}

