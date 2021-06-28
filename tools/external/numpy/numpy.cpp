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

#include <external/numpy/numpy.h>

using namespace tvb;

int searchsorted(const Vectord &a, double v, const std::string &side) {
    auto it = side.compare("left") == 0 ?
              std::lower_bound(a.begin(), a.end(), v)
                                        :
              std::upper_bound(a.begin(), a.end(), v);
    return it - a.begin();
}

Vectori searchsorted(const Vectord& a,
                              const Vectord& v,
                              const std::string& side) {
    Vectori result(v.size());
    std::transform(v.begin(), v.end(), result.begin(),
                   [&a, &side](double d) { return searchsorted(a, d, side); });
    return result;
}

double polyval(const std::vector<double>& p, double x) {
    double r = 0.0;
    unsigned N = p.size();
    for (unsigned n = 0; n < N; ++n)
        r += p[n] * pow(x, N - n - 1);
    return r;
}



Vectord tril_indices(const Matrixd& m, int N, int k) {
    int R = k < 0 ? -k : 0;
    int C = k >= 0 ? k : 0;
    Vectord y(N*N);
    int yi = 0;
    int CL = 1;
    for (int r = R; r < N; ++r) {
        for (int c = C; c < C+CL; ++c)
            y(yi++) = m(r, c);
        CL++;
    }
    return y(Eigen::seqN(0, yi));
}

tvb::Matrixd corrcoef(const tvb::Matrixd& x, const tvb::Matrixd& y, bool rowvar) {
    Matrixd c = cov(x, y, rowvar);
    Vectord d = c.matrix().diagonal().array();
    Vectord stddev = d.sqrt();
    c = c.colwise() / stddev;
    c = c.rowwise() / stddev.transpose();
    c = c.min(1.0).max(-1.0);
    return c;
}

tvb::Matrixd cov(const tvb::Matrixd& m, const tvb::Matrixd& Y, bool rowvar,
                 int ddof, bool bias) {
    Matrixd X = m;
    Matrixd y = Y;
    if (!rowvar)
        X.transposeInPlace();
    if (X.size() == 0)
        return Matrixd();
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
    Matrixd X_T = X.transpose();
    Matrixd c(X.matrix() * X_T.matrix());
    c = c * 1.0/fact;
    return c;
}

std::pair<Vectord, Vectord> average(const Matrixd& x, bool row, const Vectord& weights) {
    if (weights.size() == 0) {
        Vectord avg;
        if (row)
            avg = x.colwise().mean();
        else
            avg = x.rowwise().mean();
        double scl = double(x.size()) / avg.size();
        return { avg, Vectord::Constant(avg.size(), scl) };
    } else {
        throw("Average with weights not implemented!");
    }
}

