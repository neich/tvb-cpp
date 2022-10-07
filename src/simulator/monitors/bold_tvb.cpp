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

#include "simulator/monitors/bold_tvb.h"

using namespace tvb;


void BoldTVB::compute_hrf() {
    // Compute the hemodynamic response function
    this->m_stock_sample_rate = 0.25; // ms    # NOTE: An integral multiple of dt

    double magic_number = this->m_hrf_length; // *0.8      # truncates G, volterra kernel, once ~zero
    // Length of history needed for convolution in steps @ _stock_sample_rate
    int required_history_length = this->m_stock_sample_rate * magic_number; // 3840 for tau_s = 0.8
    this->m_stock_steps = ceil(required_history_length);
    double stock_time_max = magic_number / 1000.0; // [s]
    double stock_time_step = stock_time_max / this->m_stock_steps; // [s]
    // this->m_stock_time = TArray1d::LinSpaced(this->m_stock_steps, 0.0, stock_time_max);
    this->m_stock_time = TArray1d::LinSpaced(this->m_stock_steps, 0.0, stock_time_max - stock_time_step);
    // Compute the HRF kernel
    TArray1d G = this->m_hrf_kernel->evaluate("X", this->m_stock_time);
    m_hemodynamic_response_function = G.reverse();
    // self.hemodynamic_response_function = G[numpy.newaxis, :]
    //Interim stock configuration
    m_interim_period = 1.0 / this->m_stock_sample_rate; //period in ms
    this->m_interim_istep = int(
            round(this->m_interim_period / this->m_dt)); // interim period in integration time steps
    // self.log.debug('Bold HRF shape %s, interim period & istep %d & %d',
    //   self.hemodynamic_response_function.shape, this->m_interim_period, this->m_interim_istep)
}

void BoldTVB::config(int N, float period, float dt, const std::vector<int> &voi)  {
    init(period, dt, voi);
    m_n_nodes = N;
    m_interim_stock.clear();
    m_stock.clear();
    this->compute_hrf();
    TArray2d init_stock(m_n_nodes, m_vars_of_interest.size());
    init_stock.setZero();
    m_interim_stock.resize(m_interim_istep, init_stock);
    m_stock.resize(m_stock_steps, init_stock);
}


void BoldTVB::sample(int step, const State &state) {
    // Update the interim-stock at every step
    this->m_interim_stock[index_circ(step, m_interim_istep, -1)] = state(Eigen::all, this->m_vars_of_interest);
    // At stock's period update it with the temporal average of interim-stock


    update(step);

}

void BoldTVB::update(int step) {
    if (step % m_interim_istep == 0) {
        TArray2d avg_state(m_n_nodes, m_vars_of_interest.size());
        avg_state.setZero();
        for(auto const& is: m_interim_stock)
            avg_state += is;
        TArray2d avg_intermin_stock = avg_state / m_interim_stock.size();
        m_stock[index_circ(step / m_interim_istep, m_stock_steps, -1)] = avg_intermin_stock;
    }
    if (step % m_istep == 0) {
        int n_voi = m_vars_of_interest.size();
        int shift = (int(step / m_interim_istep) % m_stock_steps) - 1;
        TArray2d hrf = circshift(m_hemodynamic_response_function, shift);
        auto *fov = dynamic_cast<FirstOrderVolterra *>(m_hrf_kernel.get());
        State bold(m_n_nodes, n_voi);
        bold.setZero();
        if (fov != nullptr) {
            double k1_V0 = m_hrf_kernel->getVariableValue("k_1") * m_hrf_kernel->getVariableValue("V_0");
//            std::cout << "k_1 " << m_hrf_kernel->getVariableValue("k_1") << ", V_0 "
//                      << m_hrf_kernel->getVariableValue("V_0") << ", k1_V0 " << k1_V0 << std::endl;
            for (int n = 0; n < m_n_nodes; ++n)
                for (int iv = 0; iv < n_voi; ++iv) {
                    double dot = 0.0;
                    for (unsigned i = 0; i < m_stock.size(); ++i)
                        dot += m_stock[i](n, iv) * hrf(i, 0);
                    bold.col(iv)[n] = (dot - 1.0) * k1_V0;
                }
        } else {
            for (int n = 0; n < m_n_nodes; ++n)
                for (int v = 0; v < n_voi; ++v) {
                    double dot = 0.0;
                    for (unsigned i = 0; i < m_stock.size(); ++i)
                        dot += m_stock[i].col(v)[n] * hrf(i, 0);
                    bold.col(v)[n] = dot;
                }

        }
        m_records.emplace_back(Record{step * m_dt, bold});
    }
}

void BoldTVB::from_records(const std::vector<Record> &from, std::vector<Record> &to) {
    int step = 1;
    for (auto &r: from) {
        m_interim_stock[index_circ(step, m_interim_istep, -1)] = r.record;
        update(step);
    }
}
