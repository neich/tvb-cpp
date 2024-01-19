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

#ifndef TVB_CPP_BOLDTVB_H
#define TVB_CPP_BOLDTVB_H

#include <Eigen/Core>

#include "BOLDModel.h"
#include <tvb-cpp/simulator/model.h>
#include <tvb-cpp/datatypes/equation.h>

namespace tvb {

    /**
        Computes the BOLD signal from a time-series using a haemodynamic response function (HRF)

        References**:

        .. [B_1997] Buxton, R. and Frank, L., *A Model for the Coupling between
            Cerebral Blood Flow and Oxygen Metabolism During Neural Stimulation*,
            17:64-72, 1997.

        .. [Fr_2000] Friston, K., Mechelli, A., Turner, R., and Price, C., *Nonlinear
            Responses in fMRI: The Balloon Model, Volterra Kernels, and Other
            Hemodynamics*, NeuroImage, 12, 466 - 477, 2000.

        .. [Bo_1996] Geoffrey M. Boynton, Stephen A. Engel, Gary H. Glover and David
            J. Heeger (1996). Linear Systems Analysis of Functional Magnetic Resonance
            Imaging in Human V1. J Neurosci 16: 4207-4221

        .. [Po_2000] Alex Polonsky, Randolph Blake, Jochen Braun and David J. Heeger
            (2000). Neuronal activity in human primary visual cortex correlates with
            perception during binocular rivalry. Nature Neuroscience 3: 1153-1159

        .. [Gl_1999] Glover, G. *Deconvolution of Impulse Response in Event-Related BOLD fMRI*.
            NeuroImage 9, 416-429, 1999.

     */
    class BoldTVB : public BOLDModel {

        tvb::Float m_sample_rate = 0.25;
        std::unique_ptr<HRFKernelEquation> m_hrf_kernel;

        mutable tvb::Float m_hrf_length;
        mutable std::vector<TArray1d> m_stock;
        mutable int m_interim_step;
        mutable int m_istep;
        mutable int m_stock_steps;
        mutable TArray1d m_stock_time;
        mutable TArray2d m_hemodynamic_response_function;
        mutable int m_n_nodes;
        mutable tvb::Float m_dt;

        void compute_hrf() const;

        void update(int step, const State &state) const;

        void config() const;

    public:
        explicit BoldTVB(tvb::Float tr) : BOLDModel(tr), m_hrf_kernel(new FirstOrderVolterra()) {
        }

        ADD_GETTERS_AND_SETTERS_SCALAR(m_tr)

        [[nodiscard]] std::pair<TArray1d, TArray2d> compute_bold(const TArray2d &ts, tvb::Float ts_dt) const override;
    };
}


#endif //TVB_CPP_BOLDTVB_H
