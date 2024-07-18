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

#include "bold_Stephan2008.h"

#include <tvb-cpp/tools/npz_tools.h>

using namespace tvb;

TArray2d_uptr BoldStephan2008::compute_bold(const TArray2d &signal, Float dt) const {
    TArray2d signal_tmp = tvb::npy2Matrixd("signal.npy")->transpose();

    auto n_min = int(m_t_min / dt);
    auto itau = 1.0 / m_tau;
    auto ialpha = 1.0 / m_alpha;

    auto n_t = signal_tmp.cols();
    auto n_rois = signal_tmp.rows();

    TArray2d s = TArray2d::Zero(n_rois, n_t);
    TArray2d f = TArray2d::Zero(n_rois, n_t);
    TArray2d ftilde = TArray2d::Zero(n_rois, n_t);
    TArray2d v = TArray2d::Zero(n_rois, n_t);
    TArray2d vtilde = TArray2d::Zero(n_rois, n_t);
    TArray2d q = TArray2d::Zero(n_rois, n_t);
    TArray2d qtilde = TArray2d::Zero(n_rois, n_t);
    s.col(0) = 1.0;
    f.col(0) = 1.0;
    v.col(0) = 1.0;
    q.col(0) = 1.0;
    
    ftilde.col(0) = 0.0;
    vtilde.col(0) = 0.0;
    qtilde.col(0) = 0.0;

    auto dtt = dt / 1000.0; // All constants are in seconds white dt is expressed in ms
    for (int n = 0; n < n_t-1; ++n) {
        s.col(n + 1) = s.col(n) + dtt * (signal_tmp.col(n) - m_kappa * s.col(n) - m_gamma * (f.col(n) - 1.0));
        f.col(n) = f.col(n).cwiseMax(1.0);
        ftilde.col(n + 1) = ftilde.col(n) + dtt * (s.col(n) / f.col(n));
        TArray2d fv = v.col(n).pow(ialpha);

        TArray2d tmp1 = (f.col(n) - fv);
        TArray2d tmp2 = (m_tau * v.col(n));
        TArray2d tmp3 = vtilde.col(n);

        vtilde.col(n + 1) = tmp3 + dtt * (tmp1 / tmp2);
        q.col(n) = q.col(n).cwiseMax(0.01);
        TArray2d ff = (1.0 - Eigen::pow(1.0 - m_Eo, (1.0 / f.col(n)))) / m_Eo;
        qtilde.col(n + 1) = qtilde.col(n) + dtt * ((f.col(n) * ff - fv * q.col(n) / v.col(n)) / (m_tau * q.col(n)));

        f.col(n + 1) = ftilde.col(n + 1).exp();
        v.col(n + 1) = vtilde.col(n + 1).exp();
        q.col(n + 1) = qtilde.col(n + 1).exp();
        auto f_nan = tvb::isfinite(f.col(n+1));
        if (f_nan.has_value())
            printff("Hola");
        auto v_nan = tvb::isfinite(v.col(n+1));
        if (v_nan.has_value())
            printff("Hola");
        auto q_nan = tvb::isfinite(q.col(n+1));
        if (q_nan.has_value())
            printff("Hola");
    }
    tvb::Float k1 = 4.3*m_theta0*m_Eo*m_TE;
    tvb::Float k2 = m_epsilon*m_r0*m_Eo*m_TE;
    tvb::Float k3 = 1.0 - m_epsilon;
    TArray2d vv = v(Eigen::all, Eigen::seq(n_min, n_t-1));

    TArray2d qq = q(Eigen::all, Eigen::seq(n_min, n_t-1));
    vv.cwiseMin(1e-8);
    TArray2d b = m_vo * (k1 * (1 - qq) + k2 * (1 - qq / vv) + k3 * (1 - vv));
    auto step = int(1000.0 * m_tr / dt);
    auto b_size = int(b.cols());
    auto bold_size = b_size / step;
    TArray2d_uptr bold = std::make_unique<TArray2d>(n_rois, bold_size);
    int index = 0;
    for (int i = step-1; i < b_size; i += step)
        bold->col(index) = b.col(i);
    return bold;
}
