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

#include "ph_fcd.h"

#include <tvb-cpp/tools/observers/phase_interaction_matrix.h>

TArray2d PhFCD::from_fMRI(const TArray2d& signal) const {  // Compute the FCD of an input BOLD signal
    int N = signal.rows();
    int Tmax = signal.cols();
    auto npattmax = Tmax - (2*m_discard_offset-1);
    auto size_kk3 = int((npattmax - 3) * (npattmax - 2) / 2);  // The int() is not needed because N*(N-1) is always even, but "it will produce an error in the future"...

    // Isubdiag = tril_indices_column(N, k=-1)  // Indices of triangular lower part of matrix
    auto signal_filt = signal;
    if (m_apply_filters)
        signal_filt = m_filter.apply(signal);

    auto phIntMatr = phase_matrix(signal_filt, m_discard_offset);  // Compute the Phase-Interaction Matrix

    TArray1d phfcd = TArray1d::Zero(size_kk3);
    if (!tvb::isnan(signal_filt)) {
        TArrayRM2d phIntMatr_upTri = TArray2d::Zero(npattmax,
                                                  int(N * (N - 1) / 2)); // The int() is not needed, but... (see above)
        for (unsigned t = 0; t < npattmax; ++t)
            phIntMatr_upTri.row(t) = tril_values(phIntMatr[t], N, -1, true);
        int kk3 = 0;
        for (unsigned t = 0; t < npattmax - 2; ++t) {
            TVector p1 = phIntMatr_upTri(Eigen::seq(t, t + 2), Eigen::all).colwise().sum();
            auto p1_norm = p1.norm();
            for (unsigned t2 = t + 1; t2 < npattmax - 2; ++t2) {
                TVector p2 = phIntMatr_upTri.row(t2) + phIntMatr_upTri.row(t2+1) + phIntMatr_upTri.row(t2+2);;
                phfcd[kk3] = p1.dot(p2) / (p1_norm * p2.norm());
                kk3 += 1;
            }
        }
    }
    else
        phfcd = TArray1d::Constant(size_kk3, NAN);

    return phfcd;
}