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

#ifndef TVB_CPP_BOLDSTEPHAN2007_H
#define TVB_CPP_BOLDSTEPHAN2007_H

#include <Eigen/Core>

#include <simulator/monitor.h>
#include <datatypes/equation.h>
#include <eigen/circshift.h>

namespace tvb {

    class BoldStephan2007: public Monitor {
        double m_t_min = 20; // seconds

        double taus = 0.65; //  # 0.8;    # time unit (s)  --> kappa in the paper
        double tauf = 0.41; //  # 0.4;    # time unit (s)  --> gamma in the paper
        double tauo = 0.98; //  # 1;      # mean transit time (s)  --> tau in the paper
        double alpha = 0.32; // # 0.2;    # a stiffness exponent   --> alpha in the paper

        double Eo = 0.4; //  # This value is from Obata et al. (2004)
        double TE = 0.04; //  # --> TE, from Stephan et al. 2007
        double vo = 0.04; //  # ???
        double r0 = 25; //  # (s)^-1 --> r0, from Stephan et al. 2007
        double theta0 = 40.3; //  # (s)^-1

        double itaus;
        double itauf;
        double itauo;
        double ialpha;

        double k1, k2, k3;

        int m_n_nodes;
        std::vector<tvb::TArray2d> m_buffer;

        void init() {
            itaus = 1. / taus;
            itauf = 1. / tauf;
            itauo = 1. / tauo;
            ialpha = 1. / alpha;

            k1 = 4.3*theta0*Eo*TE;
            k2 = r0*Eo*TE; //  # Shouldn't it be epsilon*r0*Eo*TE ???
            k3 = 1; //  # Shouldn't it be 1-epsilon ???

            std::fill_n(m_buffer.begin(), m_vars_of_interest.size(), TArray2d::Zero(m_n_nodes, m_istep));
        }

    public:
        BoldStephan2007(int N, float period, float dt, const std::vector<int> &voi): Monitor(period, dt, voi) {
            this->config(period, N, dt, voi);
        }

        void config(float period, int N, float dt, const std::vector<int> &voi) {
            m_n_nodes = N;
            Monitor::init(period, dt, voi);
            init();
        }

        void sample(int step, const State &state) override;

    };
}


#endif //TVB_CPP_BOLD_H
