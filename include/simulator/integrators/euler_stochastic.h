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

#ifndef TVB_CPP_EULER_STOCHASTIC_H
#define TVB_CPP_EULER_STOCHASTIC_H

#include <simulator/integrator.h>
#include <simulator/noise.h>

namespace tvb {

    class EulerStochastic : public Integrator {
        Noise* m_noise;
    public:
        EulerStochastic(Noise* noise): m_noise(noise) {}

        State scheme(const State &state,
                     System &dfun,
                     const TArray2d &coupling,
                     const TArray1d &local_coupling,
                     const TArray1d &stimulus) override;

    };

}

#endif //TVB_CPP_EULER_STOCHASTIC_H
