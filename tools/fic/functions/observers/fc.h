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

#include <definitions.h>

#include <external/scipy/stats/stats.h>
#include <external/numpy/numpy.h>

#include <fic/functions/filter.h>



class FunctionalConnectivity {

protected:
    bool m_apply_filters = false;
    const Filter& m_filter;

public:

    explicit FunctionalConnectivity(bool applyFilters = false, const Filter& filter = Filter()):
    m_apply_filters(applyFilters),
    m_filter(filter)
    {}

    [[nodiscard]] virtual tvb::Matrixd from_fMRI(const tvb::Matrixd& signal) const = 0;

    virtual void init(int numSubjects, int N) = 0;

    virtual void accumulate(const tvb::Matrixd& signal, int nsub) = 0;

    [[nodiscard]] virtual double distance(const tvb::Matrixd& fcd1, const tvb::Matrixd& fcd2) const = 0;

    [[nodiscard]] virtual tvb::Matrixd postprocess() const = 0;

    // Static funcions

    static double fc_similarity(const tvb::Matrixd& fcd1, const tvb::Matrixd& fcd2) {
        Eigen::Index N = fcd1.rows();
        return FunctionalConnectivity::pearson_r(tril_indices(fcd1, N, -1), tril_indices(fcd2, N, -1));
    }

    static double pearson_r(const tvb::Vectord& x, const tvb::Vectord& y)  {
        // """Compute Pearson correlation coefficient between two arrays."""
        // Compute correlation matrix
        double corr_mat = corrcoef(x, y);
        return corr_mat;
    }

    static double KolmogorovSmirnovStatistic(const tvb::Vectord& FCD1, const tvb::Vectord& FCD2) {  // FCD similarity
        auto [d, pvalue] = ks_2samp(FCD1, FCD2);
        return d;
    }
};

class FunctionalConnectivityStandard : public FunctionalConnectivity {

protected:
    std::vector<tvb::Matrixd> m_buffer;

public:

    explicit FunctionalConnectivityStandard(bool applyFilters = false, const Filter& filter = Filter()):
            FunctionalConnectivity(applyFilters, filter)
    {}

    [[nodiscard]] virtual tvb::Matrixd from_fMRI(const tvb::Matrixd& signal) const override {
        tvb::Matrixd signal_filtered;
        if (m_apply_filters) {
            signal_filtered = m_filter.apply(signal);
            signal_filtered.transposeInPlace();
        } else {
            signal_filtered = signal.transpose();
        }
        auto cc = corrcoef(signal_filtered, tvb::Matrixd(), false);
        return cc;
    }

    virtual void init(int num_subjects, int N) override {
        m_buffer.resize(num_subjects);
        std::fill_n(m_buffer.begin(), num_subjects, tvb::Matrixd::Zero(N, N));
    }

    void accumulate(const tvb::Matrixd& signal, int nsub) override {
        m_buffer[nsub] = signal;
    }

    [[nodiscard]] double distance(const tvb::Matrixd& fcd1, const tvb::Matrixd& fcd2) const override {
        return FunctionalConnectivity::fc_similarity(fcd1, fcd2);
    }

    [[nodiscard]] tvb::Matrixd postprocess() const override {
        tvb::Matrixd result = tvb::Matrixd::Zero(m_buffer[0].rows(), m_buffer[0].cols());
        for (const tvb::Matrixd& m : m_buffer)
            result += m;
        result /= (double)m_buffer.size();
        return result;
    }
};

#endif //TVB_CPP_FC_H
