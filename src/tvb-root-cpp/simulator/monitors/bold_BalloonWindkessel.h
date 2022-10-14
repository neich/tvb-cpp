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

#include <tvb-root-cpp/simulator/monitor.h>
#include <tvb-root-cpp/datatypes/equation.h>
#include <tvb-root-cpp/eigen/circshift.h>
#include <tvb-root-cpp/simulator/simulate.h>

namespace tvb {

    class BoldBalloonWindkessel : public Monitor {

        int m_n_nodes;
        float m_t_subs;
        int m_istep_subs;
        float m_model_dt;

        float m_rho = 0.34;
        float m_alpha = 0.32;
        float m_tau = 0.98;
        float m_y = 1.0 / 0.41;
        float m_kappa = 1.0 / 0.65;
        float m_V_0 = 0.02;
        float m_k1;
        float m_k2 = 2.0;
        float m_k3;
        float m_ialpha;
        float m_itau;
        float m_oneminrho;

        TArray2d m_x0;
        TArray2d m_x1;
        TArray2d m_x2;
        TArray2d m_x3;

    public:
        BoldBalloonWindkessel(int N, float t_subs, float period, float dt, const std::vector<int> &voi): Monitor(period, dt, voi) {
            m_t_subs = t_subs;

            this->config(N, period, dt, voi);
        }

        void config(int N, float period, float dt, const std::vector<int> &voi);

        void sample(int step, const State &state) override;

    };
}


#endif //TVB_CPP_BOLDBALLONWINDKESSEL_H
