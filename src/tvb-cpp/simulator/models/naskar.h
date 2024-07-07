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

#ifndef TVB_CPP_NASKAR_H
#define TVB_CPP_NASKAR_H

#include <unsupported/Eigen/MatrixFunctions>

#include <tvb-cpp/definitions.h>
#include <tvb-cpp/simulator/model.h>

namespace tvb {

    class Naskar : public Model {
    public:
        
        TArray1d a_e;
        TArray1d a_i;
        TArray1d b_e;
        TArray1d b_i;
        TArray1d d_e;
        TArray1d d_i;
        TArray1d W_e;
        TArray1d W_i;
        TArray1d J_N;
        TArray1d gamma;
        TArray1d I0;
        TArray1d rho;
        TArray1d w;
        TArray1d t_glu;
        TArray1d t_gaba;
        TArray1d alpha_e;
        TArray1d alpha_i;
        TArray1d B_e;
        TArray1d B_i;
        TArray1d I_external;
        TArray1d M_e;
        TArray1d M_i;

    public:
        Naskar(int n_nodes) : Model(n_nodes) {
            m_cvars = { 0 };
            m_state_vars = { "S_e", "S_i", "J"};
            m_n_vars = m_state_vars.size();
            this->configure();
        }

        virtual std::vector<std::string> get_param_list() const {
            return {"a_e", "a_i",
                    "b_e", "b_i",
                    "d_e", "d_i",
                    "W_e", "W_i",
                    "J_N", "gamma", "I0", "rho", "w", "t_glu", "t_gaba", "alpha_e", "alpha_i", "B_e", "B_i",
                    "I_external", "M_e", "M_i"};
        }
        
        void set_param(const std::string& param, const TArray1d& value) override {
            ADD_SETTER_VALUE(a_e, a_i,
                             b_e, b_i,
                             d_e, d_i,
                             W_e, W_i,
                             J_N, gamma, I0, rho, w, t_glu, t_gaba, alpha_e, alpha_i, B_e, B_i,
                             I_external, M_e, M_i)
        }

        void set_param(const std::string& param, Float value) override {
            ADD_SETTER_FILL(a_e, a_i,
                            b_e, b_i,
                            d_e, d_i,
                            W_e, W_i,
                            J_N, gamma, I0, rho, w, t_glu, t_gaba, alpha_e, alpha_i, B_e, B_i,
                            I_external, M_e, M_i)
        }

        const TArray1d& get_param_value(const std::string& param) const override {
            ADD_GETTER(a_e, a_i,
                       b_e, b_i,
                       d_e, d_i,
                       W_e, W_i,
                       J_N, gamma, I0, rho, w, t_glu, t_gaba, alpha_e, alpha_i, B_e, B_i,
                       I_external, M_e, M_i)
        }

        void configure() override {

            a_e.resize(m_n_nodes);
            a_i.resize(m_n_nodes);
            b_e.resize(m_n_nodes);
            b_i.resize(m_n_nodes);
            d_e.resize(m_n_nodes);
            d_i.resize(m_n_nodes);
            W_e.resize(m_n_nodes);
            W_i.resize(m_n_nodes);
            J_N.resize(m_n_nodes);
            gamma.resize(m_n_nodes);
            I0.resize(m_n_nodes);
            rho.resize(m_n_nodes);
            w.resize(m_n_nodes);
            t_glu.resize(m_n_nodes);
            t_gaba.resize(m_n_nodes);
            alpha_e.resize(m_n_nodes);
            alpha_i.resize(m_n_nodes);
            B_e.resize(m_n_nodes);
            B_i.resize(m_n_nodes);
            I_external.resize(m_n_nodes);
            M_e.resize(m_n_nodes);
            M_i.resize(m_n_nodes);

            a_e.fill(310.0);
            a_i.fill(615.0);
            b_e.fill(125.0);
            b_i.fill(177.0);
            d_e.fill(0.16);
            d_i.fill(0.087);
            W_e.fill(1.0);
            W_i.fill(0.7);
            J_N.fill(0.15);
            gamma.fill(1.0);
            I0.fill(0.382);
            rho.fill(3.0);
            w.fill(1.4);
            t_glu.fill(7.46);
            t_gaba.fill(1.82);
            alpha_e.fill(0.072);
            alpha_i.fill(0.53);
            B_e.fill(0.0066);
            B_i.fill(0.18);
            I_external.fill(0.0);
            M_e.fill(1.0);
            M_i.fill(1.0);

        }

        [[nodiscard]] State initial() const override {
            State result(m_n_nodes, m_n_vars);
            TArray1d init_state(m_n_vars);
            init_state << 0.0, 0.0;
            for (int i = 0; i < m_n_nodes; ++i)
                result.row(i) = init_state;
            return result;
        }

        void initial(State& state) const override {
            state.resize(m_n_nodes, m_n_vars);
            TArray1d init_state(m_n_vars);
            init_state << 0.0, 0.0;
            for (int i = 0; i < m_n_nodes; ++i)
                state.row(i) = init_state;
        }

        State operator()(const State &x,
                const TArray2d &coupling,
                const TArray1d &local_coupling) const override;

    };

}

#endif //TVB_CPP_NASKAR_H
