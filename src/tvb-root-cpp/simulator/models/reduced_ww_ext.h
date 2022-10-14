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

#ifndef TVB_CPP_REDUCED_WW_EXT_H
#define TVB_CPP_REDUCED_WW_EXT_H

#include <unsupported/Eigen/MatrixFunctions>

#include <tvb-root-cpp/definitions.h>
#include <tvb-root-cpp/simulator/model.h>

namespace tvb {

    class ReducedWongWangExcInh : public Model {
    public:
        TArray1d d_i;
//    = NArray(
//            label=":external:`d_i`",
//    default=numpy.array([0.087, ]),
//    domain=Range(lo=0.01, hi=0.2, step=0.001),
//            doc="""[s]. Inhibitory population input scaling parameter chosen to fit numerical solutions.""")
//
        TArray1d a_e;
//            label=":external:`a_e`",
//    default=numpy.array([310., ]),
//    domain=Range(lo=0., hi=500., step=1.),
//            doc="[n/C]. Excitatory population input gain parameter, chosen to fit numerical solutions.")

        TArray1d b_e;
//        = NArray(
//                label = ":external:`b_e`",
//        default=numpy.array([125., ]),
//        domain = Range(lo = 0., hi = 200., step = 1.),
//                doc = "[Hz]. Excitatory population input shift parameter chosen to fit numerical solutions."
//        )

        TArray1d d_e;
//        = NArray(
//                label = ":external:`d_e`",
//        default=numpy.array([0.160, ]),
//        domain = Range(lo = 0.0, hi = 0.2, step = 0.001),
//                doc = """[s]. Excitatory population input scaling parameter chosen to fit numerical solutions."""
//        )

        TArray1d gamma_e;
//        = NArray(
//                label = r
//        ":external:`\gamma_e`",
//        default=numpy.array([0.000641, ]),
//        domain = Range(lo = 0.0, hi = 0.001, step = 0.00001),
//                doc = """Excitatory population kinetic parameter"""
//        )

        TArray1d tau_e;
//        = NArray(
//                label = r
//        ":external:`\tau_e`",
//        default=numpy.array([100., ]),
//        domain = Range(lo = 50., hi = 150., step = 1.),
//                doc = """[ms]. Excitatory population NMDA decay time constant."""
//        )

        TArray1d w_p;
//        = NArray(
//                label = r
//        ":external:`w_p`",
//        default=numpy.array([1.4, ]),
//        domain = Range(lo = 0.0, hi = 2.0, step = 0.01),
//                doc = """Excitatory population recurrence weight"""
//        )

        TArray1d J_N;
//        = NArray(
//                label = r
//        ":external:`J_{N}`",
//        default=numpy.array([0.15, ]),
//        domain = Range(lo = 0.001, hi = 0.5, step = 0.001),
//                doc = """[nA] NMDA current"""
//        )

        TArray1d W_e;
//        = NArray(
//                label = r
//        ":external:`W_e`",
//        default=numpy.array([1.0, ]),
//        domain = Range(lo = 0.0, hi = 2.0, step = 0.01),
//                doc = """Excitatory population external input scaling weight"""
//        )

        TArray1d a_i;
//        = NArray(
//                label = ":external:`a_i`",
//        default=numpy.array([615., ]),
//        domain = Range(lo = 0., hi = 1000., step = 1.),
//                doc = "[n/C]. Inhibitory population input gain parameter, chosen to fit numerical solutions."
//        )

        TArray1d b_i;
//        = NArray(
//                label = ":external:`b_i`",
//        default=numpy.array([177.0, ]),
//        domain = Range(lo = 0.0, hi = 200.0, step = 1.0),
//                doc = "[Hz]. Inhibitory population input shift parameter chosen to fit numerical solutions."
//        )

        TArray1d gamma_i;
//        = NArray(
//                label = r
//        ":external:`\gamma_i`",
//        default=numpy.array([0.001, ]),
//        domain = Range(lo = 0.0, hi = 0.002, step = 0.0001),
//                doc = """Inhibitory population kinetic parameter"""
//        )

        TArray1d tau_i;
//        = NArray(
//                label = r
//        ":external:`\tau_i`",
//        default=numpy.array([10., ]),
//        domain = Range(lo = 5., hi = 150., step = 1.0),
//                doc = """[ms]. Inhibitory population NMDA decay time constant."""
//        )

        TArray1d J_i;
//        = NArray(
//                label = r
//        ":external:`J_{i}`",
//        default=numpy.array([1, ]),
//        domain = Range(lo = 0.001, hi = 2.0, step = 0.001),
//                doc = """[nA] Local inhibitory current"""
//        )

        TArray1d W_i;
//        = NArray(
//                label = r
//        ":external:`W_i`",
//        default=numpy.array([0.7, ]),
//        domain = Range(lo = 0.0, hi = 1.0, step = 0.01),
//                doc = """Inhibitory population external input scaling weight"""
//        )

        TArray1d I_o;
//        = NArray(
//                label = ":external:`I_{o}`",
//        default=numpy.array([0.382, ]),
//        domain = Range(lo = 0.0, hi = 1.0, step = 0.001),
//                doc = """[nA]. Effective external input"""
//        )

        TArray1d G;
//        = NArray(
//                label = ":external:`G`",
//        default=numpy.array([2.0, ]),
//        domain = Range(lo = 0.0, hi = 10.0, step = 0.01),
//                doc = """Global coupling scaling"""
//        )

        TArray1d lambda;
//        = NArray(
//                label = ":external:`\lambda`",
//        default=numpy.array([0.0, ]),
//        domain = Range(lo = 0.0, hi = 1.0, step = 0.01),
//                doc = """Inhibitory global coupling scaling"""
//        )

    public:
        ReducedWongWangExcInh(int n_nodes) : Model(n_nodes) {
            this->configure(n_nodes);
            m_cvars = { 0 };
            m_state_vars = { "S_e", "S_i", "H_e", "I_e" };
            m_n_vars = m_state_vars.size();
        }

        void set_param_fill(const std::string& param, Float value) override {
            ADD_SETTER_FILL(d_i, a_e, b_e, d_e, gamma_e, tau_e, w_p, J_N, W_e, a_i, b_i, gamma_i, tau_i, J_i, W_i, I_o, G, lambda)
        }

        void set_param_value(const std::string& param, const TArray1d& value) override {
            ADD_SETTER_VALUE(d_i, a_e, b_e, d_e, gamma_e, tau_e, w_p, J_N, W_e, a_i, b_i, gamma_i, tau_i, J_i, W_i, I_o, G, lambda)
        }

        const TArray1d& get_param(const std::string& param) const override {
            ADD_GETTER(d_i, a_e, b_e, d_e, gamma_e, tau_e, w_p, J_N, W_e, a_i, b_i, gamma_i, tau_i, J_i, W_i, I_o, G, lambda)
        }



        void configure(int n_nodes) {
            d_i.resize(n_nodes);
            a_e.resize(n_nodes);
            b_e.resize(n_nodes);
            d_e.resize(n_nodes);
            gamma_e.resize(n_nodes);
            tau_e.resize(n_nodes);
            w_p.resize(n_nodes);
            J_N.resize(n_nodes);
            W_e.resize(n_nodes);
            a_i.resize(n_nodes);
            b_i.resize(n_nodes);
            gamma_i.resize(n_nodes);
            tau_i.resize(n_nodes);
            J_i.resize(n_nodes);
            W_i.resize(n_nodes);
            I_o.resize(n_nodes);
            G.resize(n_nodes);
            lambda.resize(n_nodes);


            d_i.fill(0.087);
            a_e.fill(310.0);
            b_e.fill(125.0);
            d_e.fill(0.16);
            gamma_e.fill(0.000641);
            tau_e.fill(100.0);
            w_p.fill(1.4);
            J_N.fill(0.15);
            W_e.fill(1.0);
            a_i.fill(615.0);
            b_i.fill(177.0);
            gamma_i.fill(0.001);
            tau_i.fill(10.0);
            J_i.fill(1.0);
            W_i.fill(0.7);
            I_o.fill(0.382);
            G.fill(1.0);
            lambda.fill(0.0);
        }

        [[nodiscard]] State initial() const override {
            State result(m_n_nodes, m_n_vars);
            TArray1d init_state(m_n_vars);
            init_state << 0.0, 0.0, 0.0, 0.0;
            for (int i = 0; i < m_n_nodes; ++i)
                result.row(i) = init_state;
            return result;
        }

        void initial(State& state) const override {
            state.resize(m_n_nodes, m_n_vars);
            TArray1d init_state(m_n_vars);
            init_state << 0.001, 0.001, 0.0, 0.0;
            for (int i = 0; i < m_n_nodes; ++i)
                state.row(i) = init_state;
        }



        State operator()(const State &x,
                const TArray2d &coupling,
                const TArray1d &local_coupling) override;

    };

}

#endif //TVB_CPP_REDUCED_WW_EXT_H
