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

#define _USE_MATH_DEFINES
#include <cmath>

#include "phase_interaction_matrix.h"

#include <tvb-cpp/tools/scipy/signal/signaltools.h>
#include <tvb-cpp/tools/numpy/numpy.h>


std::vector<TArray2d> tvb::phase_matrix(const TArray2d& signal, unsigned discard_offset) {
    int N = signal.rows();
    int Tmax = signal.cols();
    auto npattmax = Tmax - (2*discard_offset-1);
    auto PhIntMatr = std::vector<TArray2d>(npattmax, TArray2d::Zero(N, N));

    if (!tvb::isnan(signal)) {
        // Data structures we are going to need...
        TArray2d phases = TArray2d::Zero(N, Tmax);
        TArray2d dFC = TArray2d::Zero(N, N);
        // PhIntMatr = np.zeros((npattmax, int(N * (N - 1) / 2)))  # The int() is not needed, but... (see above)
        // syncdata = np.zeros(npattmax)

        for (unsigned n = 0; n < N; ++n) {
            TArray1dc Xanalytic = hilbert_signal(signal.row(n) - signal.row(n).mean());
            phases.row(n) = Xanalytic.array().arg();
        }

        // Isubdiag = tril_indices_column(N, k=-1)  # Indices of triangular lower part of matrix
        for (unsigned t = discard_offset; t < Tmax - discard_offset + 1; ++t) {
            // kudata = np.sum(np.cos(phases[:, t - 1]) + 1j * np.sin(phases[:, t - 1])) / N
            // syncdata[t - 10] = abs(kudata)
            for (unsigned i = 0; i < N; ++i)
                for (unsigned j = 0; j < i+1; ++j) {
                    // print(f'processing {t}: ({i}, {j})')
                    tvb::Float da = fabs(phases(i, index_circ(t-1, Tmax)) - phases(j, index_circ(t-1, Tmax)));
                    auto dcos = cos(da > M_PI ? 2.0 * M_PI - da : da);
                    // da = fmod(da + 180.0, 360.0) - 180.0;
                    dFC(i, j) = dcos;
                    dFC(j, i) = dcos;
                }
            PhIntMatr[index_circ(t - discard_offset, Tmax)] = dFC;

        }
    }
    else {
        throw std::runtime_error("tvb::phase_matrix produced a result that contains nan!");
    }

    return PhIntMatr;
}