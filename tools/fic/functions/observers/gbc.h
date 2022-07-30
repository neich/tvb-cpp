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

#ifndef TVB_CPP_GBC_H
#define TVB_CPP_GBC_H

#include <definitions.h>

#include <external/scipy/stats/stats.h>
#include <external/numpy/numpy.h>

#include <fic/functions/filter.h>
#include <fic/functions/observers/fc.h>


class GBC_FC : public FunctionalConnectivityStandard {

public:

    explicit GBC_FC(bool applyFilters = false, const Filter& filter = Filter()):
            FunctionalConnectivityStandard(applyFilters, filter)
    {}

    [[nodiscard]] double distance(const tvb::TArray2d& fcd1, const tvb::TArray2d& fcd2) const override {
        return FunctionalConnectivity::pearson_r(fcd1.row(0), fcd2.row(0));
    }

    [[nodiscard]] tvb::TArray2d postprocess() const override {
        tvb::TArray2d fc_emp = FunctionalConnectivityStandard::postprocess();
        int N = fc_emp.rows();
        tvb::TArray2d onesd = TArray2d::Zero(N, N);
        onesd.matrix().diagonal() = TArray1d::Ones(N);
        tvb::TArray2d result = fc_emp * onesd;
        result = fc_emp - result;
        result = result.colwise().mean();
        return result;
    }
};

#endif //TVB_CPP_GBC_H
