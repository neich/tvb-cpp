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

#include "bold_Stephan2007.h"

using namespace tvb;

TArray2d_uptr BoldStephan2007::compute_bold(const TArray2d &ts_bold, Float ts_dt) const {
    TArray2d ts = ts_bold.transpose(); // Transpose to get memory contiguous vector operations

    int n_t = ts.cols(); // number of time samples

    auto N = ts.rows(); // number of ROIs

    TArray2d x0 = TArray2d::Zero(N, n_t);
    TArray2d x1 = TArray2d::Ones(N, n_t);
    TArray2d x2 = TArray2d::Ones(N, n_t);
    TArray2d x3 = TArray2d::Ones(N, n_t);
    tvb::Float dt = .001 * ts_dt; // The next formulas use dt in seconds and ts_dt is given in ms
    for (unsigned n = 0; n < n_t - 1; ++n) {
        //  # Shouldn't it be (0.5 r[n] + 3) instead of r[n] ??? also, shouldn't it be m_taus and m_tauf instead of itaus and itauf???
        x0.col(n + 1) = x0.col(n) + dt * (ts.col(n) - itaus * x0.col(n) - itauf * (x1.col(n) - 1.0));
// # Equation (10) for f in Stephan et al. 2007
        x1.col(n + 1) = x1.col(n) + dt * x0.col(n);
//# Equation (8) for v and q in Stephan et al. 2007
        x2.col(n + 1) = x2.col(n) + dt * itauo * (x1.col(n) - x2.col(n).pow(ialpha));
        x3.col(n + 1) = x3.col(n) + dt * itauo * (x1.col(n) * (1.0 - pow((1 - m_eo), (1 / x1.col(n)))) / m_eo -
                                                    x2.col(n).pow(ialpha) * x3.col(n) / x2.col(n));
    }
    TArray2d v = x2; // (Eigen::all, Eigen::seq(n_min, m_istep - 1));
    TArray2d q = x3; // (Eigen::all, Eigen::seq(n_min, m_istep - 1));
    TArray2d b = m_vo * (k1 * (Float(1.0) - q) + k2 * (Float(1.0) - q / v) + k3 * (Float(1.0) - v));

    return std::make_unique<TArray2d>(b.transpose());
}
