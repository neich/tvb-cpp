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

#include "bold_filters.h"

using namespace tvb;


TArray2d BandPassFilter::apply(const TArray2d &boldSignal) const {
    int N = boldSignal.rows();
    // int Tmax = boldSignal.cols();
    double fnq = 1. / (2. * m_TR); //             # Nyquist frequency
    TArray1d Wn(2);
    Wn[0] = m_flp / fnq;
    Wn[1] = m_fhi / fnq; // # butterworth bandpass non-dimensional frequency
    auto[bfilt, afilt] = iirfilter_ba(k, Wn); //   # construct the filter
    TArray2d signal_filt = TArray2d::Zero(boldSignal.rows(), boldSignal.cols());
    for (int seed = 0; seed < N; ++seed) {
        TArray1d ts = detrend_linear(boldSignal.row(seed));
        ts -= ts.mean();
        double std_dev = std::sqrt((ts - ts.mean()).square().sum() / (ts.size() - 1));
        for (unsigned i = 0; i < ts.size(); ++i) {
            if (ts[i] > 3. * std_dev) ts[i] = 3. * std_dev;    // # Remove strong artefacts
            if (ts[i] < -3. * std_dev) ts[i] = -3. * std_dev;    // # Remove strong artefacts
        }
        signal_filt.row(seed) = filtfilt_pad(bfilt, vc2vd(afilt), ts,
                                             3 * (std::max(bfilt.size(), afilt.size()) - 1)); //  # Band pass filter. padlen modified to get the same result as in Matlab
    }
    return signal_filt;

}
