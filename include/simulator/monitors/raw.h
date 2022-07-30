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

#ifndef TVB_CPP_RAW_H
#define TVB_CPP_RAW_H

#include <simulator/monitor.h>

namespace tvb {

    class Raw: public Monitor {
    public:

        void config_for_sim(const SimConfig &sim_config, const TArray1di &voi) override {
            if (this->period() != sim_config.dt())
                this->m_dt = sim_config.dt();
            Monitor::config_for_sim(sim_config, voi);
            this->m_istep = 1;
        }

        virtual MSample sample(int step, const State &state) override {
            return MSample(step * this->m_dt, state);
        }
    };
}

#endif //TVB_CPP_RAW_H
