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

#include <tvb-cpp/tools/scipy/stats/stats.h>

#include "fc.h"


class PhFCD : public FunctionalConnectivity {
protected:
    tvb::TArray2d_uptr m_buffer;
    unsigned m_discard_offset = 0;

public:
    explicit PhFCD(unsigned discard_offset = 0, bool applyFilters = false, Filter_uptr filter = nullptr) :
            FunctionalConnectivity(applyFilters, std::move(filter)),
            m_discard_offset(discard_offset) {}

    tvb::TArray2d_uptr from_fMRI(const tvb::TArray2d &signal) const override;

    void init(int numSubjects, int N) override {
        m_buffer = std::make_unique<tvb::TArray2d>();
    }

    void accumulate(const tvb::TArray2d &signal, int nsub = 0) override {
        assert(("Signal for phFCD has to be a vector (1 column)!", signal.cols() == 1));
        tvb::TArray2d_uptr new_buffer = tvb::TArray2d_uptr(new tvb::TArray2d(m_buffer->rows() + signal.rows(), 1));
        if (m_buffer->size() > 0)
            *new_buffer << *m_buffer, signal.col(0);
        m_buffer = std::move(new_buffer);
    }

    TArray2d_uptr postprocess() override {
        auto uptr = std::move(m_buffer);
        return uptr;
    }

    double distance(const tvb::TArray2d &fcd1, const tvb::TArray2d &fcd2) const override {
        return FunctionalConnectivity::KolmogorovSmirnovStatistic(fcd1.col(0), fcd2.col(0));
    }

};