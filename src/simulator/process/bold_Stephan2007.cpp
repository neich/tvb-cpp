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

#include "simulator/monitors/bold_Stephan2007.h"

using namespace tvb;


//MSample BoldStephan2007::sample(int step, const State &state) {
//
//    return MSample(-1, {});
//}

TArray2d BoldStephan2007::apply(const TArray2d& signal) {
    // # Initial conditions

    TArray2d x0 = TArray2d::Zero(signal.rows(), signal.cols());
    TArray2d x1 = TArray2d::Zero(signal.rows(), signal.cols());
    TArray2d x2 = TArray2d::Zero(signal.rows(), signal.cols());
    TArray2d x3 = TArray2d::Zero(signal.rows(), signal.cols());
    x0.col(0) = 0.0;
    x1.col(0) = 1.0;
    x2.col(0) = 1.0;
    x3.col(0) = 1.0;

    int n_t = int(m_T/m_dt);

    int n_min = int(round(m_t_min / m_dt));

    for (unsigned n = 0; n < n_t-1; ++n) {
        //  # Shouldn't it be (0.5 r[n] + 3) instead of r[n] ??? also, shouldn't it be taus and tauf instead of itaus and itauf???
        x0.col(n + 1) = x0.col(n) + m_dt * (signal.col(n) - itaus * x0.col(n) - itauf * (x1.col(n) - 1.0));
// # Equation (10) for f in Stephan et al. 2007
        x1.col(n + 1) = x1.col(n) + m_dt * x0.col(n);
//# Equation (8) for v and q in Stephan et al. 2007
        x2.col(n + 1) = x2.col(n) + m_dt * itauo * (x1.col(n) - x2.col(n).pow(ialpha));
        x3.col(n + 1) = x3.col(n) + m_dt * itauo * (x1.col(n) * (1.0 - pow((1 - Eo), (1/x1.col(n))))/Eo - x2.col(n).pow(ialpha) * x3.col(n)/x2.col(n));
    }
    TArray2d v = x2(Eigen::all, Eigen::seq(n_min, n_t - 1));
    TArray2d q = x3(Eigen::all, Eigen::seq(n_min, n_t - 1));
    TArray2d b = vo * (k1 * (1 - q) + k2 * (1 - q / v) + k3 * (1 - v));


    return b;
}
