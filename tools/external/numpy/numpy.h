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

#ifndef TVB_CPP_NUMPY_H
#define TVB_CPP_NUMPY_H

#include <string>
#include <vector>

#include <definitions.h>

int searchsorted(const tvb::Vectord &a, double v, const std::string &side);

tvb::Vectori searchsorted(const tvb::Vectord &a, const tvb::Vectord &v, const std::string &side);

double polyval(const std::vector<double> &p, double x);


template<typename _Scalar>
Eigen::Array<_Scalar, Eigen::Dynamic,1> arange(_Scalar start, _Scalar end, _Scalar step = 1) {
    assert(step != 0);
    _Scalar d = abs(end - start);
    int N = int(d / step);
    if (step >= 1.0 && int(d) % int(step) != 0) N++;
    Eigen::Array<_Scalar, Eigen::Dynamic,1> result(N);
    for (int i = 0; i < N; ++i, start+=step)
        result[i] = start;
    return result;
}

inline
tvb::Vectord cumsum(const tvb::Vectord& a) {
    tvb::Vectord result(a.size());
    double accum = 0.0;
    for (int i = 0; i != a.size(); ++i) {
        result[i] = accum + a[i];
        accum += a[i];
    }

    return result;

}

tvb::Matrixd corrcoef(const tvb::Matrixd& x,
                      const tvb::Matrixd& y=tvb::Matrixd(),
                      bool rowvar=true);

tvb::Matrixd cov(const tvb::Matrixd& x,
                 const tvb::Matrixd& y=tvb::Matrixd(),
                 bool rowvar=true,
                 int ddof=-1,
                 bool bias=false);

std::pair<tvb::Vectord, tvb::Vectord> average(const tvb::Matrixd& x,
                                        bool row=true,
                                        const tvb::Vectord& weights={});

tvb::Vectord tril_indices(const tvb::Matrixd& m, int N, int k=0);

#endif //TVB_CPP_NUMPY_H
