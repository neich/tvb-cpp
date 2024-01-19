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


class SW_FC : public FunctionalConnectivity {
protected:
    int m_windowSize = 30;
    int m_windowStep = 3;
    tvb::TArray1d m_buffer;

public:
    explicit SW_FC(bool applyFilters = false, const Filter &filter = Filter()) : FunctionalConnectivity(applyFilters, filter) {}

    SW_FC(int windowSize, int windowStep, bool applyFilters = false, const Filter &filter = Filter()) :
            FunctionalConnectivity(applyFilters, filter),
            m_windowSize(windowSize), m_windowStep(windowStep) {}

    tvb::TArray2d from_fMRI(const tvb::TArray2d &signal) const override;

    void init(int numSubjects, int N) override {
        m_buffer = tvb::TArray1d();
    }

    void accumulate(const tvb::TArray2d& signal, int nsub = 0)  override {
        tvb::TArray1d new_buffer = tvb::TArray1d(m_buffer.size() + signal.rows());
        new_buffer << m_buffer, signal.col(0);
        m_buffer = new_buffer;
    }

    tvb::TArray2d postprocess() const override {
        return m_buffer;
    }

    double distance(const tvb::TArray2d &fcd1, const tvb::TArray2d &fcd2) const override {
        return FunctionalConnectivity::KolmogorovSmirnovStatistic(fcd1.col(0), fcd2.col(0));
    }

};