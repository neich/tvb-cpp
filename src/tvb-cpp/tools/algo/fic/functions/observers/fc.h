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

#ifndef TVB_CPP_FC_H
#define TVB_CPP_FC_H

#include <tvb-cpp/definitions.h>

#include <tvb-cpp/tools/algo/external/scipy/stats/stats.h>
#include <tvb-cpp/tools/algo/external/numpy/numpy.h>

#include <tvb-cpp/tools/algo/fic/functions/filter.h>



class FunctionalConnectivity {

protected:
    bool m_apply_filters = false;
    const Filter& m_filter;

public:

    explicit FunctionalConnectivity(bool applyFilters = false, const Filter& filter = Filter()):
    m_apply_filters(applyFilters),
    m_filter(filter)
    {}

    [[nodiscard]] virtual tvb::TArray2d from_fMRI(const tvb::TArray2d& signal) const = 0;

    virtual void init(int numSubjects, int N) = 0;

    virtual void accumulate(const tvb::TArray2d& signal, int nsub = 0) = 0;

    [[nodiscard]] virtual double distance(const tvb::TArray2d& fcd1, const tvb::TArray2d& fcd2) const = 0;

    [[nodiscard]] virtual tvb::TArray2d postprocess() const = 0;

    // Static funcions

    static double fc_similarity(const tvb::TArray2d& fcd1, const tvb::TArray2d& fcd2) {
        Eigen::Index N = fcd1.rows();
        return FunctionalConnectivity::pearson_r(tril_indices(fcd1, N, -1), tril_indices(fcd2, N, -1));
    }

    static double pearson_r(const tvb::TArray1d& x, const tvb::TArray1d& y)  {
        // """Compute Pearson correlation coefficient between two arrays."""
        // Compute correlation matrix
        double corr_mat = corrcoef(x, y);
        return corr_mat;
    }

    static double KolmogorovSmirnovStatistic(const tvb::TArray1d& FCD1, const tvb::TArray1d& FCD2) {  // FCD similarity
        auto [d, pvalue] = ks_2samp(FCD1, FCD2);
        return d;
    }
};

class FunctionalConnectivityStandard : public FunctionalConnectivity {

protected:
    std::vector<tvb::TArray2d> m_buffer;

public:

    explicit FunctionalConnectivityStandard(bool applyFilters = false, const Filter& filter = Filter()):
            FunctionalConnectivity(applyFilters, filter)
    {}

    [[nodiscard]] virtual tvb::TArray2d from_fMRI(const tvb::TArray2d& signal) const override {
        tvb::TArray2d signal_filtered;
        if (m_apply_filters) {
            signal_filtered = m_filter.apply(signal);
            signal_filtered.transposeInPlace();
        } else {
            signal_filtered = signal.transpose();
        }
        auto cc = corrcoef(signal_filtered, tvb::TArray2d(), false);
        return cc;
    }

    virtual void init(int num_subjects, int N) override {
        m_buffer.resize(num_subjects);
        std::fill_n(m_buffer.begin(), num_subjects, tvb::TArray2d::Zero(N, N));
    }

    void accumulate(const tvb::TArray2d& signal, int nsub) override {
        m_buffer[nsub] = signal;
    }

    [[nodiscard]] double distance(const tvb::TArray2d& fcd1, const tvb::TArray2d& fcd2) const override {
        return FunctionalConnectivity::fc_similarity(fcd1, fcd2);
    }

    [[nodiscard]] tvb::TArray2d postprocess() const override {
        tvb::TArray2d result = tvb::TArray2d::Zero(m_buffer[0].rows(), m_buffer[0].cols());
        for (const tvb::TArray2d& m : m_buffer)
            result += m;
        result /= (double)m_buffer.size();
        return result;
    }
};

#endif //TVB_CPP_FC_H
