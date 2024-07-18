//
// Created by natxm on 08/12/2023.
//

#ifndef TVB_CPP_BOLDMODEL_H
#define TVB_CPP_BOLDMODEL_H

#include <tvb-cpp/param_object.h>

namespace tvb {

    class BOLDModel : public ParamObject {
    protected:
        tvb::Float m_tr;

    public:
        explicit BOLDModel(tvb::Float tr = 1.0) : m_tr(tr) {}

        /**

        * @param ts time series as a TArray2d(n_time_samples, n_rois)
        * @param ts_dt time series time step
        * @return pair of TArray1d(n_time_samples_bold, n_rois) and TArray2d(n_time_samples, n_rois)
        *
        * n_time_samples_bold is usually <= n_time_samples
        */
        [[nodiscard]] virtual TArray2d_uptr compute_bold(const tvb::TArray2d &ts, tvb::Float ts_dt) const = 0;

    };

}
#endif //TVB_CPP_BOLDMODEL_H
