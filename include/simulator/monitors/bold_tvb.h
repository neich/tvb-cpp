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

#ifndef TVB_CPP_BOLDTVB_H
#define TVB_CPP_BOLDTVB_H

#include <Eigen/Core>

#include <simulator/monitor.h>
#include <datatypes/equation.h>
#include <eigen/circshift.h>
#include "simulator/simulate.h"

namespace tvb {

    class BoldTVB : public Monitor {

        int m_n_nodes;
        double m_hrf_length = 20000.0;
        std::unique_ptr<HRFKernelEquation> m_hrf_kernel;

        double m_interim_period;
        int m_interim_istep;
//
        std::vector<State> m_interim_stock;
        std::vector<State> m_stock;
        int m_stock_steps;
        TArray1d m_stock_time;
        double m_stock_sample_rate = 0.02;
        TArray2d m_hemodynamic_response_function;

        void compute_hrf();

    public:
        BoldTVB(int N, float period, float dt, const std::vector<int> &voi): Monitor(period, dt, voi), m_hrf_kernel(new FirstOrderVolterra()) {
            this->config(N, period, dt, voi);
        }

        void config(int N, float period, float dt, const std::vector<int> &voi);

        void sample(int step, const State &state) override;

    };
}


#endif //TVB_CPP_BOLD_H
