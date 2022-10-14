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

#ifndef TVB_CPP_MONTBRIO_H
#define TVB_CPP_MONTBRIO_H

#include <unsupported/Eigen/MatrixFunctions>

#include <tvb-root-cpp/definitions.h>
#include <tvb-root-cpp/simulator/model.h>

namespace tvb {

    class Montbrio : public Model {
    public:

        TArray1d tau_e;
//        = NArray(
//                label = ":external:`\tau`",
//        default=numpy.array([10., ]),
//        domain = Range(lo = 1., hi = 50., step = 1.),
//                doc = """[ms]. Excitatory population decay time constant."""
//        )

        TArray1d tau_i;
//        = NArray(
//                label = ":external:`\tau`",
//        default=numpy.array([10., ]),
//        domain = Range(lo = 1., hi = 50., step = 1.),
//                doc = """[ms]. Inhibitory population decay time constant."""
//        )

        TArray1d G;
//        = NArray(
//                label = ":external:`G`",
//        default=numpy.array([2.0, ]),
//        domain = Range(lo = 0.0, hi = 10.0, step = 0.01),
//                doc = """Global coupling scaling"""
//        )


        TArray1d delta_e;
        TArray1d delta_i;
        TArray1d eta_e;
        TArray1d eta_i;
        TArray1d I_e;
        TArray1d I_i;
        TArray1d I_ext;
        Float J_ee;
        Float J_ei;
        Float J_ie;
        Float J_ii;
        TArray1d J_e;
        TArray1d J_i;

    protected:
        Float m_dt;
        int m_num_steps;

    public:
        Montbrio(int n_nodes, Float t_start, Float t_end, Float dt) : Model(n_nodes) {
            m_dt = dt;
            m_num_steps = floor((t_end - t_start) / m_dt);
            this->configure(n_nodes);
            m_cvars = {4}; // coupling variable S_ee
            m_state_vars = {"r_e", "r_i", "u_e", "u_i", "S_ee", "S_ei", "S_ie", "S_ii"};
            m_n_vars = m_state_vars.size();
        }

        void set_param(const std::string &param, Float value) {
            ADD_SETTER_FILL(tau_e, tau_i, G, delta_e, delta_i, eta_e, eta_i, I_e, I_i, I_ext, J_e, J_i)
        }

        void configure(int n_nodes) {
            tau_e.resize(n_nodes);
            tau_i.resize(n_nodes);
            G.resize(n_nodes);
            delta_e.resize(n_nodes);
            delta_i.resize(n_nodes);
            eta_e.resize(n_nodes);
            eta_i.resize(n_nodes);
            I_e.resize(n_nodes);
            I_i.resize(n_nodes);
            I_ext.resize(n_nodes);
            J_e.resize(n_nodes);
            J_i.resize(n_nodes);


            tau_e.fill(10.0);
            tau_i.fill(10.0);
            J_i.fill(0.0);
            G.fill(2.5);
            delta_e.fill(1.0);
            delta_i.fill(1.0);
            eta_e.fill(1.0);
            eta_i.fill(1.0);
            I_e.fill(0.0);
            I_i.fill(0.0);
            I_ext.fill(0.0);
            J_e.fill(1.0);
            J_i.fill(0.0);
            J_ee = 1.0;
            J_ei = 0.0;
            J_ie = 0.0;
            J_ii = 0.0;

            init_dependant();
        }

        void init_dependant() override {
            // TODO tau_av = floor(Float(1e-3) * tau / m_dt + Float(1e-6));
//            m_I = TArray1d::Zero(m_num_steps);
//            m_utrace = TArray2d::Zero(m_n_nodes, m_num_steps+1);
//            m_r = TArray2d::Zero(m_n_nodes, m_num_steps+1);
        }


        [[nodiscard]] State initial() const override {
            State result(m_n_nodes, m_n_vars);
            TArray1d init_state(m_n_vars);
            init_state << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
            for (int i = 0; i < m_n_nodes; ++i)
                result.row(i) = init_state;
            return result;
        }

        void initial(State &state) const override {
            state.resize(m_n_nodes, m_n_vars);
            TArray1d init_state(m_n_vars);
            init_state << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
            for (int i = 0; i < m_n_nodes; ++i)
                state.row(i) = init_state;
        }


        State operator()(const State &x,
                         const TArray2d &coupling,
                         const TArray1d &local_coupling) override;

    };

}

#endif //TVB_CPP_MONTBRIO_H
