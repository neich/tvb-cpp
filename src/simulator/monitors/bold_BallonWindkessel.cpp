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

#include "simulator/monitors/bold_BalloonWindkessel.h"

using namespace tvb;

void BoldBalloonWindkessel::config(int N, float period, float dt, const std::vector<int> &voi) {
    init(period, dt, voi);
    m_k1 = 7 * m_rho;
    m_k3 = 2 * m_rho - 0.2;
    m_ialpha = 1.0 / m_alpha;
    m_itau = 1.0 / m_tau;
    m_oneminrho = 1.0 - m_rho;
    m_n_nodes = N;
    m_model_dt = m_dt / 100;
    m_istep_subs = lround(m_t_subs / m_dt);

    m_x0 = TArray2d::Constant(m_n_nodes, m_vars_of_interest.size(), 0.0);
    m_x1 = TArray2d::Constant(m_n_nodes, m_vars_of_interest.size(), 1.0);
    m_x2 = TArray2d::Constant(m_n_nodes, m_vars_of_interest.size(), 1.0);
    m_x3 = TArray2d::Constant(m_n_nodes, m_vars_of_interest.size(), 1.0);
}


void BoldBalloonWindkessel::sample(int step, const State &state) {
    // Update the interim-stock at every step
    if (step % m_istep_subs == 0) {
        for (int voi = 0; voi < m_vars_of_interest.size(); ++voi) {
            m_x0.col(voi) = m_x0.col(voi) + m_model_dt * (state.col(m_vars_of_interest[voi]) - m_kappa * m_x0.col(voi) - m_y * (m_x1.col(voi) - 1.0));
            auto ftmp = m_x1.col(voi) + m_model_dt * m_x0.col(voi);
            m_x2.col(voi) = m_x2.col(voi) + m_model_dt * m_itau * (m_x1.col(voi) - m_x2.col(voi).pow(m_ialpha));
            m_x3.col(voi) = m_x3.col(voi) + m_model_dt * m_itau * (m_x1.col(voi) * (1.0 - Eigen::pow(m_oneminrho, (1.0 / m_x1.col(voi)))) / m_rho -
                    m_x2.col(voi).pow(m_ialpha) * m_x3.col(voi) / m_x2.col(voi));
            m_x1.col(voi) = ftmp;
        }
    }

    if (step % m_istep == 0) {
        TArray2d bold = 100.0 / m_rho * m_V_0 *
                        (m_k1 * (1.0 - m_x3) + m_k2 * (1.0 - m_x3 / m_x2) +
                         m_k3 * (1 - m_x2));
        m_records.push_back(Monitor::Record{step * m_dt, bold});
    }

}
