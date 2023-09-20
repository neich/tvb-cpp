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

#include "numpy.h"

using namespace tvb;

int searchsorted(const TArray1d &a, double v, const std::string &side) {
    auto it = side.compare("left") == 0 ?
              std::lower_bound(a.begin(), a.end(), v)
                                        :
              std::upper_bound(a.begin(), a.end(), v);
    return it - a.begin();
}

TArray1di searchsorted(const TArray1d& a,
                       const TArray1d& v,
                       const std::string& side) {
    TArray1di result(v.size());
    std::transform(v.begin(), v.end(), result.begin(),
                   [&a, &side](double d) { return searchsorted(a, d, side); });
    return result;
}

Float polyval(const std::vector<Float>& p, Float x) {
    double r = 0.0;
    unsigned N = p.size();
    for (unsigned n = 0; n < N; ++n)
        r += p[n] * pow(x, N - n - 1);
    return r;
}



TArray1d tril_indices(const TArray2d& m, int N, int k) {
    int R = k < 0 ? -k : 0;
    int C = k >= 0 ? k : 0;
    TArray1d y(N * N);
    int yi = 0;
    int CL = 1;
    for (int r = R; r < N; ++r) {
        for (int c = C; c < C+CL; ++c)
            y(yi++) = m(r, c);
        CL++;
    }
    return y(Eigen::seqN(0, yi));
}

tvb::TArray2d corrcoef(const tvb::TArray2d& x, const tvb::TArray2d& y, bool rowvar) {
    TArray2d c = cov(x, y, rowvar);
    TArray1d d = c.matrix().diagonal().array();
    TArray1d stddev = d.sqrt();
    c = c.colwise() / stddev;
    c = c.rowwise() / stddev.transpose();
    c = c.min(1.0).max(-1.0);
    return c;
}

tvb::TArray2d cov(const tvb::TArray2d& m, const tvb::TArray2d& Y, bool rowvar,
                  int ddof, bool bias) {
    TArray2d X = m;
    TArray2d y = Y;
    if (!rowvar)
        X.transposeInPlace();
    if (X.size() == 0)
        return TArray2d();
    if (y.size() > 0)  {
        if (!rowvar)
            y.transposeInPlace();
        X << y;
    }
    if (ddof < 0) {
        if (bias)
            ddof = 0;
        else
            ddof = 1;
    }

    auto [avg, w_sum_v] = average(X, false);
    // double w_sum = w_sum_v[0];
    int fact = X.cols() - ddof;
    if (fact <= 0)
        fact = 0;

    X = X.colwise() - avg;
    TArray2d X_T = X.transpose();
    TArray2d c(X.matrix() * X_T.matrix());
    c = c * 1.0/fact;
    return c;
}

std::pair<TArray1d, TArray1d> average(const TArray2d& x, bool row, const TArray1d& weights) {
    if (weights.size() == 0) {
        TArray1d avg;
        if (row)
            avg = x.colwise().mean();
        else
            avg = x.rowwise().mean();
        double scl = double(x.size()) / avg.size();
        return {avg, TArray1d::Constant(avg.size(), scl) };
    } else {
        throw("Average with weights not implemented!");
    }
}

