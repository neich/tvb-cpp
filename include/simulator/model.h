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

        void append(const StateTrack &state_track) {
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

        virtual void initial(State &) const {
            throw std::runtime_error("Model initial state not implemented!");
        }

        virtual StateTrack *create_track() const {
            return new StateTrack();
        }

        const std::vector<int> &cvars() const { return m_cvars; }

        int n_vars() const { return m_n_vars; }

        const std::vector<std::string> &state_vars() const { return m_state_vars; }

        State operator()(const State &x,
                         const TArray2d &coupling,
                         const TArray1d &local_coupling) const = 0;

    };


    #define FE_0(WHAT)
    #define FE_1(WHAT, X) WHAT(X)
    #define FE_2(WHAT, X, ...) WHAT(X)FE_1(WHAT, __VA_ARGS__)
    #define FE_3(WHAT, X, ...) WHAT(X)FE_2(WHAT, __VA_ARGS__)
    #define FE_4(WHAT, X, ...) WHAT(X)FE_3(WHAT, __VA_ARGS__)
    #define FE_5(WHAT, X, ...) WHAT(X)FE_4(WHAT, __VA_ARGS__)
    #define FE_6(WHAT, X, ...) WHAT(X)FE_5(WHAT, __VA_ARGS__)
    #define FE_7(WHAT, X, ...) WHAT(X)FE_6(WHAT, __VA_ARGS__)
    #define FE_8(WHAT, X, ...) WHAT(X)FE_7(WHAT, __VA_ARGS__)
    #define FE_9(WHAT, X, ...) WHAT(X)FE_8(WHAT, __VA_ARGS__)
    #define FE_10(WHAT, X, ...) WHAT(X)FE_9(WHAT, __VA_ARGS__)
    #define FE_11(WHAT, X, ...) WHAT(X)FE_10(WHAT, __VA_ARGS__)
    #define FE_12(WHAT, X, ...) WHAT(X)FE_11(WHAT, __VA_ARGS__)
    #define FE_13(WHAT, X, ...) WHAT(X)FE_12(WHAT, __VA_ARGS__)
    #define FE_14(WHAT, X, ...) WHAT(X)FE_13(WHAT, __VA_ARGS__)
    #define FE_15(WHAT, X, ...) WHAT(X)FE_14(WHAT, __VA_ARGS__)
    #define FE_16(WHAT, X, ...) WHAT(X)FE_15(WHAT, __VA_ARGS__)
    #define FE_17(WHAT, X, ...) WHAT(X)FE_16(WHAT, __VA_ARGS__)
    #define FE_18(WHAT, X, ...) WHAT(X)FE_17(WHAT, __VA_ARGS__)
    #define FE_19(WHAT, X, ...) WHAT(X)FE_18(WHAT, __VA_ARGS__)
    #define FE_20(WHAT, X, ...) WHAT(X)FE_19(WHAT, __VA_ARGS__)
    #define FE_21(WHAT, X, ...) WHAT(X)FE_20(WHAT, __VA_ARGS__)
    #define FE_22(WHAT, X, ...) WHAT(X)FE_21(WHAT, __VA_ARGS__)
    #define FE_23(WHAT, X, ...) WHAT(X)FE_22(WHAT, __VA_ARGS__)
    #define FE_24(WHAT, X, ...) WHAT(X)FE_23(WHAT, __VA_ARGS__)
    #define FE_25(WHAT, X, ...) WHAT(X)FE_24(WHAT, __VA_ARGS__)
    #define FE_26(WHAT, X, ...) WHAT(X)FE_25(WHAT, __VA_ARGS__)
    #define FE_27(WHAT, X, ...) WHAT(X)FE_26(WHAT, __VA_ARGS__)
    #define FE_28(WHAT, X, ...) WHAT(X)FE_27(WHAT, __VA_ARGS__)
    #define FE_29(WHAT, X, ...) WHAT(X)FE_28(WHAT, __VA_ARGS__)
    #define FE_30(WHAT, X, ...) WHAT(X)FE_29(WHAT, __VA_ARGS__)

    #define GET_MACRO(_0,_1,_2,_3,_4,_5,_6,_7,_8,_9,_10,_11,_12,_13,_14,_15,_16,_17,_18,_19,_20,_21,_22,_23,_24,_25,_26,_27,_28,_29,_30,NAME,...) NAME
    #define FOR_EACH(action,...) \
      GET_MACRO(_0,__VA_ARGS__,FE_30,FE_29,FE_28,FE_27,FE_26,FE_25,FE_24,FE_23,FE_22,FE_21,FE_20,FE_19,FE_18,FE_17,FE_16,FE_15,FE_14,FE_13,FE_12,F_11,FE_10,FE_9,FE_8,FE_7,FE_6,FE_5,FE_4,FE_3,FE_2,FE_1,FE_0)(action,__VA_ARGS__)

    #define SETTER(field) if (#field == param) this->field.fill(value); \

    #define ADD_SETTER(...) FOR_EACH(SETTER, __VA_ARGS__)
}

#endif //TVB_CPP_MODEL_H
