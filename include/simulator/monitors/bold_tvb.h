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

namespace tvb {

    class BoldTVB : public Monitor {

        double m_hrf_length = 20000.0;
        std::unique_ptr<HRFKernelEquation> m_hrf_kernel;

        double m_interim_period;
        int m_interim_istep;
//
        std::vector<State> m_interim_stock;
        std::vector<State> m_stock;
        int m_stock_steps;
        Vectord m_stock_time;
        double m_stock_sample_rate = 0.02;
        Matrixd m_hemodynamic_response_function;

        void compute_hrf();

        void init();

    public:
        BoldTVB(const SimConfig &sim_config, const std::vector<int> &voi) : m_hrf_kernel(new FirstOrderVolterra()) {
            m_period = 0.5;
            this->config_for_sim(sim_config, voi);
        }
        BoldTVB(int N, double dt, const std::vector<int> &voi) : m_hrf_kernel(new FirstOrderVolterra()) {
            m_period = 0.5;
            this->config_for_sim(N, dt, voi);
        }

        void config_for_sim(const SimConfig &sim_config, const std::vector<int> &voi) override {
            Monitor::config_for_sim(sim_config, voi);
            init();
        }

        void config_for_sim(int N, double dt, const std::vector<int> &voi) override {
            Monitor::config_for_sim(N, dt, voi);
            init();
        }

        MSample sample(int step, const State &state) override;

    };
}


#endif //TVB_CPP_BOLD_H
