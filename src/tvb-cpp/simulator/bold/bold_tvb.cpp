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

#include "bold_tvb.h"
#include <tvb-cpp/eigen/circshift.h>

using namespace tvb;


void BoldTVB::compute_hrf() const {
    // Compute the hemodynamic response function

    // Compute the HRF kernel
    tvb::Float required_history_length = m_sample_rate * m_hrf_length;
    m_stock_steps = std::ceil(required_history_length);
    tvb::Float time_max = m_hrf_length / 1000.0;
    tvb::Float time_step = time_max / m_stock_steps;
    TArray1d time = tvb::arange(0.0, time_max, time_step);
    TArray1d G = this->m_hrf_kernel->evaluate("X", time);
    m_hemodynamic_response_function = G.reverse();
    // self.hemodynamic_response_function = G[numpy.newaxis, :]
    //Interim stock configuration
    int interim_period = int(1.0 / m_sample_rate);
    m_interim_step = int(round(interim_period / m_dt));
    this->m_istep = int(
            round(this->m_tr / this->m_dt)); // interim period in integration time steps
    // self.log.debug('Bold HRF shape %s, interim period & istep %d & %d',
    //   self.hemodynamic_response_function.shape, this->m_interim_period, this->m_interim_istep)
}

void BoldTVB::config() const {
    m_stock.clear();
    this->compute_hrf();
    TArray2d init_stock(m_n_nodes, 1);
    init_stock.setZero();
    m_stock.resize(m_stock_steps, init_stock);
}

void BoldTVB::update(int step, const State &state) const {
}

std::pair<TArray1d, TArray2d> BoldTVB::compute_bold(const TArray2d &ts, tvb::Float ts_dt) const {
    m_n_nodes = ts.cols(); // number of ROIsn
    int n_samples = ts.rows();
    m_dt = ts_dt;
    m_hrf_length = m_tr * 10.0;
    this->config();

    int n_bold = n_samples / m_istep;
    TArray2d data(n_bold, m_n_nodes);
    TArray1d t_samples(n_bold);
    int index = 0;
    for (int step = 1; step <= n_samples; ++step) {
        if (step % m_interim_step == 0) {
            int start = step-m_interim_step;
            int finish = step;
            TArray1d a_sample = TArray1d::Zero(m_n_nodes);
            for (int i = start; i < finish; ++i) {
                const TArray1d& sample = ts.row(i);
                a_sample += sample;
            }
            a_sample /= tvb::Float(m_interim_step);
            this->m_stock[(step / m_interim_step % m_stock_steps) - 1] = a_sample;
        }
        if (step % m_istep == 0) {
            int shift = int(step / m_interim_step % m_stock_steps) - 1;
            TArray2d hrf = circshift(m_hemodynamic_response_function, shift);
            auto *fov = dynamic_cast<FirstOrderVolterra *>(m_hrf_kernel.get());
            TArray1d bold = TArray1d::Zero(m_n_nodes);
            if (fov != nullptr) {
                double k1_V0 = m_hrf_kernel->getVariableValue("k_1") * m_hrf_kernel->getVariableValue("V_0");
//            std::cout << "k_1 " << m_hrf_kernel->getVariableValue("k_1") << ", V_0 "
//                      << m_hrf_kernel->getVariableValue("V_0") << ", k1_V0 " << k1_V0 << std::endl;
                for (int n = 0; n < m_n_nodes; ++n) {
                    double dot = 0.0;
                    for (unsigned i = 0; i < m_stock.size(); ++i)
                        dot += m_stock[i](n) * hrf(i, 0);
                    bold[n] = (dot - 1.0) * k1_V0;
                }
            } else {
                for (int n = 0; n < m_n_nodes; ++n) {
                        double dot = 0.0;
                        for (unsigned i = 0; i < m_stock.size(); ++i)
                            dot += m_stock[i][n] * hrf(i, 0);
                        bold[n] = dot;
                    }

            }
            t_samples[index] = step * m_dt;
            data.row(index) = bold;
            index++;
        }
    }
    return {t_samples, data};
}
