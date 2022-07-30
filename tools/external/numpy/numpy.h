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

int searchsorted(const tvb::TArray1d &a, double v, const std::string &side);

tvb::TArray1di searchsorted(const tvb::TArray1d &a, const tvb::TArray1d &v, const std::string &side);

tvb::Float polyval(const std::vector<tvb::Float> &p, tvb::Float x);


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
tvb::TArray1d cumsum(const tvb::TArray1d& a) {
    tvb::TArray1d result(a.size());
    double accum = 0.0;
    for (int i = 0; i != a.size(); ++i) {
        result[i] = accum + a[i];
        accum += a[i];
    }

    return result;

}

tvb::TArray2d corrcoef(const tvb::TArray2d& x,
                       const tvb::TArray2d& y=tvb::TArray2d(),
                       bool rowvar=true);

tvb::TArray2d cov(const tvb::TArray2d& x,
                  const tvb::TArray2d& y=tvb::TArray2d(),
                  bool rowvar=true,
                  int ddof=-1,
                  bool bias=false);

std::pair<tvb::TArray1d, tvb::TArray1d> average(const tvb::TArray2d& x,
                                                bool row=true,
                                                const tvb::TArray1d& weights={});

tvb::TArray1d tril_indices(const tvb::TArray2d& m, int N, int k=0);

#endif //TVB_CPP_NUMPY_H
