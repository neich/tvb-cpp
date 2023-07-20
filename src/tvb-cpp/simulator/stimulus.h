//
// Created by imartin on 26-Oct-22.
//

#ifndef TVB_ROOT_CPP_PROJECTS_STIMULUS_H
#define TVB_ROOT_CPP_PROJECTS_STIMULUS_H


#include <tvb-cpp/simulator/model.h>

namespace tvb {

    class Stimulus {
    protected:
        Float m_t_start;
        Float m_t_end;
        Float m_dt;
        std::vector<int> m_vois;

    public:
        Stimulus() = default;
        Stimulus(const std::vector<int>& vois): m_vois(vois) {}

        void configure(Float t_start, Float t_end, Float dt) {
            m_t_start = t_start;
            m_t_end = t_end;
            m_dt = dt;
        }

        virtual ~Stimulus() = default;

        virtual State update(int step, const State &stimulus) = 0;
        virtual State initial(State state) = 0;
    };

    class NullStimulus : public Stimulus {
    public:
        NullStimulus() = default;

        State update(int step, const State &stimulus) override {
            return stimulus;
        }
        State initial(State state) override {
            return State::Zero(state.rows(), state.cols());
        }

    };
}
#endif //TVB_ROOT_CPP_PROJECTS_STIMULUS_H
