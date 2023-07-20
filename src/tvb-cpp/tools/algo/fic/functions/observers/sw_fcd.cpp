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

#include <tvb-cpp/simulator/monitors/bold_tvb.h>

#include <tvb-cpp/tools/algo/external/scipy/stats/stats.h>
#include <tvb-cpp/tools/algo/external/numpy/numpy.h>
#include <tvb-cpp/tools/algo/fic/functions/observers/sw_fcd.h>


int calc_length(int start, int end, int step) {
    // This fails for a negative step e.g., range(10, 0, -1).
    // From https://stackoverflow.com/questions/31839032/python-how-to-calculate-the-length-of-a-range-without-creating-the-range
    return 1 + (end - start - 1) / step;
}

TArray2d SW_FC::from_fMRI(const TArray2d& signal) const {  // Compute the FCD of an input BOLD signal
    int N = signal.rows();
    int Tmax = signal.cols();
    TArray2d signal_filtered;
    if (m_apply_filters)
        signal_filtered = m_filter.apply(signal);  // Filters seem to be always applied...
    else
        signal_filtered = signal;

//    Isubdiag = np.tril_indices(N, k=-1)  // Indices of triangular lower part of matrix

    // For each pair of sliding windows calculate the FC at t and t2 and
    // compute the correlation between the two.
    int lastWindow = Tmax - m_windowSize;  // 190 = 220 - 30
    int N_windows = calc_length(0, lastWindow, m_windowStep);  // N_windows = len(np.arange(0, lastWindow, windowStep))
    TArray1d cotsampling = TArray1d::Zero((int(N_windows * (N_windows - 1) / 2)));
    int kk = 0;
    int ii2 = 0;
    for (int t = 0; t < lastWindow; t+=m_windowStep) {
        int jj2 = 0;
        TArray2d sfilt = signal_filtered(Eigen::all, Eigen::seqN(t, m_windowSize + 1)).transpose();
        TArray2d cc = corrcoef(sfilt, TArray2d(), false);  // Pearson correlation coefficients
        for (int t2 = 0; t2 < lastWindow; t2+=m_windowStep) {
            TArray2d sfilt2 = signal_filtered(Eigen::all, Eigen::seqN(t2, m_windowSize + 1)).transpose();
            TArray2d cc2 = corrcoef(sfilt2, TArray2d(), false);  // Pearson correlation coefficients
            double ca = pearson_r(tril_indices(cc, N, -1), tril_indices(cc2, N, -1));  // Correlation between both FC
            if (jj2 > ii2) {  // Only keep the upper triangular part
                cotsampling[kk++] = ca;
            }
            jj2++;
        }
        ii2++;
    }
    return cotsampling;

}