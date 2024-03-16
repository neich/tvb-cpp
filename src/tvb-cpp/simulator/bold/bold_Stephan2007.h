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

#ifndef TVB_CPP_BOLDSTEPHAN2007_H
#define TVB_CPP_BOLDSTEPHAN2007_H

#include "BOLDModel.h"

namespace tvb {

    class BoldStephan2007 : public BOLDModel {
    protected:
        tvb::Float m_t_min = 20000; // miliseconds

        tvb::Float m_taus = 0.65; //  # 0.8;    # time unit (s)  --> kappa in the paper
        tvb::Float m_tauf = 0.41; //  # 0.4;    # time unit (s)  --> gamma in the paper
        tvb::Float m_tauo = 0.98; //  # 1;      # mean transit time (s)  --> tau in the paper
        tvb::Float m_alpha = 0.32; // # 0.2;    # a stiffness exponent   --> m_alpha in the paper

        tvb::Float m_eo = 0.4; //  # This value is from Obata et al. (2004)
        tvb::Float m_te = 0.04; //  # --> m_te, from Stephan et al. 2007
        tvb::Float m_vo = 0.04; //  # ???
        tvb::Float m_r0 = 25; //  # (s)^-1 --> m_r0, from Stephan et al. 2007
        tvb::Float m_theta0 = 40.3; //  # (s)^-1

        tvb::Float itaus{};
        tvb::Float itauf{};
        tvb::Float itauo{};
        tvb::Float ialpha{};

        tvb::Float k1{}, k2{}, k3{};

        virtual void init() {
            itaus = 1. / m_taus;
            itauf = 1. / m_tauf;
            itauo = 1. / m_tauo;
            ialpha = 1. / m_alpha;

            k1 = 4.3 * m_theta0 * m_eo * m_te;
            k2 = m_r0 * m_eo * m_te; //  # Shouldn't it be epsilon*m_r0*m_eo*m_te ???
            k3 = 1; //  # Shouldn't it be 1-epsilon ???
        }

    public:
        explicit BoldStephan2007(tvb::Float tr = 1.0) : BOLDModel(tr) {
            this->init();
        }

        ADD_GETTERS_AND_SETTERS_SCALAR(m_taus, m_tauf, m_tauo, m_alpha, m_eo, m_te, m_vo, m_r0, m_theta0, m_tr)

        [[nodiscard]] std::pair<TArray1d, TArray2d> compute_bold(const TArray2d &ts, tvb::Float ts_dt) const override;
    };


    class BoldStephan2007b : public BoldStephan2007 {

    public:
        explicit BoldStephan2007b(tvb::Float tr = 1.0) : BoldStephan2007(tr) {
            m_eo = 0.4;
            m_te = 0.04;
            m_vo = 0.04;
            m_theta0 = 40.3;

            this->init_dependant();
        }

        void init_dependant() override {
            itaus = 1. / m_taus;
            itauf = 1. / m_tauf;
            itauo = 1. / m_tauo;
            ialpha = 1. / m_alpha;

            k1 = 4.3 * 40.3 * m_theta0 * m_eo * m_te;
            k2 = 25.0 * m_eo * m_te;
            k3 = 1;
        }

        ADD_GETTERS_AND_SETTERS_SCALAR(m_taus, m_tauf, m_tauo, m_alpha, m_eo, m_te, m_vo, m_theta0, m_tr)
    };
}


#endif //TVB_CPP_BOLDSTEPHAN2007_H
