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

#ifndef TVB_CPP_BOLDSTEPHAN2008_H
#define TVB_CPP_BOLDSTEPHAN2008_H

#include "BOLDModel.h"

namespace tvb {

    class BoldStephan2008 : public BOLDModel {
    protected:
        tvb::Float m_t_min = 20000; // miliseconds

        tvb::Float m_kappa = 0.65;  // 0.8;    // Rate of vasodilatory signal decay, time unit (s) [Friston2003], tvb::Float m_eta = 0.64; in [Friston2019]
        tvb::Float m_gamma = 0.41;  // 0.4;    // Rate of flow-dependent elimination, time unit (s)  [Friston2003], tvb::Float m_chi = 0.32; in [Friston2019]
        tvb::Float m_tau = 0.98;  // 1;      // mean transit time (s) in [Friston2003], 1/tvb::Float m_tau = 2; in [Friston2019]
        tvb::Float m_alpha = 0.32; //0.32; % 0.2;    % Grubb's exponent (a stiffness exponent) [Friston2003] and [Friston2019]
        tvb::Float m_epsilon = 0.34; // Intravascular:extravascular ratio (should be one epsilon for each brain area.)
// [Buxton et al. 1998] used 0.4...
// Lu and Van Zijl 2005 found values that make it 1 [Stephan et al. 2007]...
// [Friston2019] initializes this as 1

        tvb::Float m_Eo = 0.4;   // region-specific resting oxygen extraction fractions. This value is from [Obata et al. 2004]
        tvb::Float m_TE = 0.04;  // Echo time (seconds), TE, from [Stephan et al. 2007].
// In [Friston2019] they use phi = Eo*TE as another variable
        tvb::Float m_vo = 0.08;  // resting venous volume, DCM CODE reads 100*0.08...
        tvb::Float m_r0 = 25;  // (s)^-1 --> slope r0 of intravascular relaxation rate R_iv as a function of oxygen
// saturation Y:  R_iv = r0*[(1-Y)-(1-Y0)]. This value of r0 from [Stephan et al. 2007]
        tvb::Float m_theta0 = 40.3;  // (s)^-1, frequency offset at the outer surface of magnetized vessels


        virtual void init() {
        }

    public:
        explicit BoldStephan2008(tvb::Float tr = 1.0) : BOLDModel(tr) {
            this->init();
        }

        ADD_GETTERS_AND_SETTERS_SCALAR(m_kappa, m_gamma, m_tau, m_alpha, m_epsilon, m_Eo, m_TE, m_vo, m_r0, m_theta0)

        [[nodiscard]] TArray2d_uptr compute_bold(const TArray2d &signal, Float ts_dt) const override;
    };


}


#endif //TVB_CPP_BOLDSTEPHAN2008_H
