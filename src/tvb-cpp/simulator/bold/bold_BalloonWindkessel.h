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

#ifndef TVB_CPP_BOLDBALLONWINDKESSEL_H
#define TVB_CPP_BOLDBALLONWINDKESSEL_H

#include <Eigen/Core>

#include <tvb-cpp/simulator/monitor.h>
#include <tvb-cpp/datatypes/equation.h>
#include <tvb-cpp/eigen/circshift.h>
#include <tvb-cpp/simulator/simulate.h>

namespace tvb {

    class BoldBalloonWindkessel : public Monitor {

        int m_n_nodes;
        tvb::Float m_t_subs;
        int m_istep_subs;
        tvb::Float m_model_dt;

        tvb::Float m_rho = 0.34;
        tvb::Float m_alpha = 0.32;
        tvb::Float m_tau = 0.98;
        tvb::Float m_y = 1.0 / 0.41;
        tvb::Float m_kappa = 1.0 / 0.65;
        tvb::Float m_V_0 = 0.02;
        tvb::Float m_k1;
        tvb::Float m_k2 = 2.0;
        tvb::Float m_k3;
        tvb::Float m_ialpha;
        tvb::Float m_itau;
        tvb::Float m_oneminrho;

        TArray2d m_x0;
        TArray2d m_x1;
        TArray2d m_x2;
        TArray2d m_x3;

    public:
        BoldBalloonWindkessel(int N, tvb::Float t_subs, tvb::Float period, tvb::Float dt, const std::vector<int> &voi): Monitor(period, dt, voi) {
            m_t_subs = t_subs;

            this->config(N, period, dt, voi);
        }

        void config(int N, tvb::Float period, tvb::Float dt, const std::vector<int> &voi);

        void sample(int step, const State &state) override;

    };
}


#endif //TVB_CPP_BOLDBALLONWINDKESSEL_H
