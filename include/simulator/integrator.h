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

#ifndef TVB_CPP_INTEGRATOR_H
#define TVB_CPP_INTEGRATOR_H

#include <memory>
#include <vector>

// #include <boost/ref.hpp>
// #include <boost/numeric/odeint.hpp>

#include <simulator/model.h>


namespace tvb {


    struct StoreStateTrack {
        std::vector<State> &m_states;
        std::vector<double> &m_times;

        StoreStateTrack(std::vector<State> &states, std::vector<double> &times) : m_states(states), m_times(times) {}

        void operator()(const State &x, double t) {
            m_states.push_back(x);
            m_times.push_back(t);
        }

    };


    class Integrator {
    protected:
        double m_start_time;
        double m_end_time;
        double m_dt;
    public:

        typedef typename std::unique_ptr<Integrator> UPtr;

        Integrator() {
            this->configure(0.0, 1000.0, 0.1);
        }

        Integrator(double start_time, double end_time, double dt) {
            this->configure(start_time, end_time, dt);
        }

        void configure(double start_time, double end_time, double dt) {
            m_start_time = start_time;
            m_end_time = end_time;
            m_dt = dt;
        }

        double dt() const { return m_dt; }

        virtual State scheme(const State &state,
                             System &sys,
                             const TArray2d &coupling,
                             const TArray1d &local_coupling,
                             const TArray1d &stimulus) = 0;
    };
}

#endif //TVB_CPP_INTEGRATOR_H
