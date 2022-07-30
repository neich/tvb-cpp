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

#ifndef TVB_CPP_MODEL_H
#define TVB_CPP_MODEL_H

#include <memory>
#include <vector>

#include <definitions.h>

namespace tvb {


    typedef TArray2d State; // (node, variable)

    class StateTrack {
    public:
        std::vector<State> m_states;
        std::vector<double> m_times;

        const std::vector<State> &states() const { return m_states; }

        std::vector<State> &states() { return m_states; }

        const std::vector<double> &times() const { return m_times; }

        std::vector<double> &times() { return m_times; }

        const State &latest() { return m_states.back(); }

        void push(const State &state, double time) {
            m_states.push_back(state);
            m_times.push_back(time);
        }

        void append(const StateTrack& state_track) {
            m_states.reserve(state_track.states().size());
            m_times.reserve(state_track.times().size());
            m_states.insert(m_states.end(), state_track.states().begin(), state_track.states().end());
            m_times.insert(m_times.end(), state_track.times().begin(), state_track.times().end());
        }

    };

    class System {
    public:

        virtual State operator()(const State &x, const TArray2d &coupling, const TArray1d &local_coupling) const = 0;
    };


    class Model : public System {
    protected:
        int m_n_nodes;
        int m_n_vars;

        std::vector<int> m_cvars;
        std::vector<std::string> m_state_vars;

    public:
        typedef typename std::unique_ptr<Model> UPtr;


        Model(int n_nodes, int n_vars) : m_n_nodes(n_nodes), m_n_vars(n_vars) {
        }

        virtual State initial() const {
            throw std::runtime_error("Model initial state not implemented!");
        }

        virtual void initial(State&) const {
            throw std::runtime_error("Model initial state not implemented!");
        }

        virtual StateTrack *create_track() const {
            return new StateTrack();
        }

        const std::vector<int>& cvars() const { return m_cvars; }
        int n_vars() const { return m_n_vars; }

        const std::vector<std::string>& state_vars() const { return  m_state_vars; }

        State operator()(const State &x,
                const TArray2d &coupling,
                const TArray1d &local_coupling) const = 0;

    };
}

#endif //TVB_CPP_MODEL_H
