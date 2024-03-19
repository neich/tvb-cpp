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

#include <tvb-cpp/definitions.h>
#include <tvb-cpp/simulator/model.h>

namespace tvb {

    class Montbrio : public Model {
    public:

        TArray1d tau_e;
        TArray1d tau_i;
        TArray1d tau_N;
        TArray1d delta_e;
        TArray1d delta_i;
        TArray1d eta_e;
        TArray1d eta_i;
        TArray1d g_e;
        TArray1d g_i;
        TArray1d g_ee;
        TArray1d g_ei;
        TArray1d g_ie;
        TArray1d g_ii;

        TArray1d a_e;
        TArray1d a_i;
        TArray1d I_e_ext;
        TArray1d I_i_ext;
        TArray1d J_e;
        TArray1d J_i;
        TArray1d J_ee;
        TArray1d J_ei;
        TArray1d J_ie;
        TArray1d J_ii;
        TArray1d J_A;
        TArray1d J_N_ee;
        TArray1d J_N_ie;
        TArray1d J_G_ei;
        TArray1d J_G_ii;
        TArray1d J;

    protected:
        int m_num_steps;

    public:
        Montbrio(int m_n_nodes) : Model(m_n_nodes) {
            m_cvars = {0}; // coupling variable r_e
            m_state_vars = {"r_e", "r_i", "u_e", "u_i", "S_ee", "S_ie"};
            m_n_vars = m_state_vars.size();
            this->configure();
        }

        std::vector<std::string> get_param_list() const override {
            return {"tau_e", "tau_i", "tau_N", "delta_e", "delta_i", "eta_e", "eta_i", "a_e", "a_i", "g_e", "g_i",
                    "I_e_ext", "I_i_ext", "J_e", "J_i", "J_A", "J_G_ei", "J_G_ii", "J_N_ee", "J_N_ie", "J"};
        }
        
        void set_param(const std::string& param, const TArray1d& value) override {
            ADD_SETTER_VALUE(tau_e, tau_i, tau_N, delta_e, delta_i, eta_e, eta_i, a_e, a_i, g_e, g_i, I_e_ext, I_i_ext,\
                             J_e, J_i, J_A, J_G_ei, J_G_ii, J_N_ee, J_N_ie, J)
            throw std::runtime_error(string_format("ParamArray %s does not exist in this model", param.c_str()));

        }

        void set_param(const std::string& param, Float value) override {
            ADD_SETTER_FILL(tau_e, tau_i, tau_N, delta_e, delta_i, eta_e, eta_i, a_e, a_i, g_e, g_i, I_e_ext, I_i_ext,\
                            J_e, J_i, J_A, J_G_ei, J_G_ii, J_N_ee, J_N_ie, J)
            throw std::runtime_error(string_format("ParamScalar %s does not exist in this model", param.c_str()));
        }

        const TArray1d& get_param_value(const std::string& param) const override {
            ADD_GETTER(tau_e, tau_i, tau_N, delta_e, delta_i, eta_e, eta_i, a_e, a_i, g_e, g_i, I_e_ext, I_i_ext, J_e,\
                       J_i, J_A, J_G_ei, J_G_ii, J_N_ee, J_N_ie, J)
            throw std::runtime_error(string_format("ParamScalar %s does not exist in this model", param.c_str()));
        }

        void configure() override {
            tau_e.resize(m_n_nodes);
            tau_i.resize(m_n_nodes);
            tau_N.resize(m_n_nodes);
            delta_e.resize(m_n_nodes);
            delta_i.resize(m_n_nodes);
            eta_e.resize(m_n_nodes);
            eta_i.resize(m_n_nodes);
            a_e.resize(m_n_nodes);
            a_i.resize(m_n_nodes);
            g_e.resize(m_n_nodes);
            g_i.resize(m_n_nodes);
            g_ee.resize(m_n_nodes);
            g_ei.resize(m_n_nodes);
            g_ie.resize(m_n_nodes);
            g_ii.resize(m_n_nodes);

            I_e_ext.resize(m_n_nodes);
            I_i_ext.resize(m_n_nodes);
            J_e.resize(m_n_nodes);
            J_i.resize(m_n_nodes);
            J_A.resize(m_n_nodes);
            J_ee.resize(m_n_nodes);
            J_ei.resize(m_n_nodes);
            J_ie.resize(m_n_nodes);
            J_ii.resize(m_n_nodes);

            J_G_ei.resize(m_n_nodes);
            J_G_ii.resize(m_n_nodes);
            J_N_ee.resize(m_n_nodes);
            J_N_ie.resize(m_n_nodes);
            J.resize(m_n_nodes);

            tau_e.fill(10.0);
            tau_i.fill(10.0);
            tau_N.fill(10.0);
            delta_e.fill(1.0);
            delta_i.fill(1.0);
            eta_e.fill(1.0);
            eta_i.fill(1.0);
            a_e.fill(0.25);
            a_i.fill(1.0);
            g_e.fill(2.5);
            g_i.fill(0);
            g_ee.fill(2.5);
            g_ei.fill(0.0);
            g_ie.fill(2.5);
            g_ii.fill(0.0);

            I_e_ext.fill(0.0);
            I_i_ext.fill(0.0);
            J_e.fill(1.0);
            J_i.fill(0.0);
            J_A.fill(1.0);
            J_ee.fill(10.0);
            J_ei.fill(10.0);
            J_ie.fill(10.0);
            J_ii.fill(10.0);
            J.fill(10.0);

            J_G_ei.fill(1.0);
            J_G_ii.fill(1.0);
            J_N_ee.fill(1.0);
            J_N_ie.fill(1.0);
        }

        void init_dependant() override {
            J_N_ee = J_ee + g_ee*a_e.log();
            J_N_ie = J_ie + g_ie*a_e.log();
            J_G_ei = J_ei + g_ei*a_i.log();
            J_G_ii = J_ii + g_ii*a_i.log();
        }

        [[nodiscard]] State initial() const override {
            State result(m_n_nodes, m_n_vars);
            TArray1d init_state = TArray1d::Zero(m_n_vars);
            for (int i = 0; i < m_n_nodes; ++i)
                result.row(i) = init_state;
            return result;
        }

        void initial(State &state) const override {
            state.resize(m_n_nodes, m_n_vars);
            TArray1d init_state = TArray1d::Zero(m_n_vars);
            for (int i = 0; i < m_n_nodes; ++i)
                state.row(i) = init_state;
        }


        State operator()(const State &x,
                         const TArray2d &coupling,
                         const TArray1d &local_coupling) const override;

    };

    /*class MontbrioSimple : public Model {
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

        TArray1d tau_N;

        TArray1d delta_e;
        TArray1d delta_i;
        TArray1d eta_e;
        TArray1d eta_i;
        TArray1d g_e;
        TArray1d g_i;
        TArray1d a_e;
        TArray1d a_i;
        TArray1d I_e_ext;
        TArray1d I_i_ext;
        TArray1d J_e;
        TArray1d J_i;
        TArray1d J_A;
        TArray1d J_G;
        TArray1d J_N;
        TArray1d J;

    protected:
        int m_num_steps;

    public:
        MontbrioSimple(int m_n_nodes) : Model(m_n_nodes) {
            m_cvars = {0}; // coupling variable r_e
            m_state_vars = {"r_e", "r_i", "u_e", "u_i", "S_e"};
            m_n_vars = m_state_vars.size();
        }

        virtual std::vector<std::string> get_param_list() const {
            return {"tau_e", "tau_i", "tau_N", "delta_e", "delta_i", "eta_e", "eta_i", "a_e", "a_i", "g_e", "g_i", "I_e_ext", "I_i_ext", "J_e", "J_i", "J_A", "J_G", "J_N", "J"};
        }

        void set_param(const std::string& param, const TArray1d& value) override {
            ADD_SETTER_VALUE(tau_e, tau_i, tau_N, delta_e, delta_i, eta_e, eta_i, a_e, a_i, g_e, g_i, I_e_ext, I_i_ext, J_e, J_i, J_A, J_G, J_N, J)
            throw std::runtime_error(string_format("ParamArray %s does not exist in this model", param.c_str()));

        }

        void set_param(const std::string& param, Float value) override {
            ADD_SETTER_FILL(tau_e, tau_i, tau_N, delta_e, delta_i, eta_e, eta_i, a_e, a_i, g_e, g_i, I_e_ext, I_i_ext, J_e, J_i, J_A, J_G, J_N, J)
            throw std::runtime_error(string_format("ParamScalar %s does not exist in this model", param.c_str()));
        }

        const TArray1d& get_param_value(const std::string& param) const override {
            ADD_GETTER(tau_e, tau_i, tau_N, delta_e, delta_i, eta_e, eta_i, a_e, a_i, g_e, g_i, I_e_ext, I_i_ext, J_e, J_i, J_A, J_G, J_N, J)
            throw std::runtime_error(string_format("ParamScalar %s does not exist in this model", param.c_str()));
        }

        void configure() override {
            tau_e.resize(m_n_nodes);
            tau_i.resize(m_n_nodes);
            tau_N.resize(m_n_nodes);
            delta_e.resize(m_n_nodes);
            delta_i.resize(m_n_nodes);
            eta_e.resize(m_n_nodes);
            eta_i.resize(m_n_nodes);
            a_e.resize(m_n_nodes);
            a_i.resize(m_n_nodes);
            g_e.resize(m_n_nodes);
            g_i.resize(m_n_nodes);
            I_e_ext.resize(m_n_nodes);
            I_i_ext.resize(m_n_nodes);
            J_e.resize(m_n_nodes);
            J_i.resize(m_n_nodes);
            J_A.resize(m_n_nodes);
            J_G.resize(m_n_nodes);
            J_N.resize(m_n_nodes);
            J.resize(m_n_nodes);

            tau_e.fill(10.0);
            tau_i.fill(10.0);
            tau_N.fill(10.0);
            delta_e.fill(1.0);
            delta_i.fill(1.0);
            eta_e.fill(1.0);
            eta_i.fill(1.0);
            a_e.fill(0.25);
            a_i.fill(1.0);
            g_e.fill(2.5);
            g_i.fill(0);
            I_e_ext.fill(0.0);
            I_i_ext.fill(0.0);
            J_e.fill(1.0);
            J_i.fill(0.0);
            J_A.fill(1.0);
            J_G.fill(1.0);
            J_N.fill(1.0);
            J.fill(10.0);
        }

        void init_dependant() override {
            J_N = J_e + g_e*a_e.log();
            J_G = J_i + g_i*a_i.log();
        }

        [[nodiscard]] State initial() const override {
            State result(m_n_nodes, m_n_vars);
            TArray1d init_state(m_n_vars);
            init_state << 0.0, 0.0, 0.0, 0.0, 0.0;
            for (int i = 0; i < m_n_nodes; ++i)
                result.row(i) = init_state;
            return result;
        }

        void initial(State &state) const override {
            state.resize(m_n_nodes, m_n_vars);
            TArray1d init_state(m_n_vars);
            init_state << 0.0, 0.0, 0.0, 0.0, 0.0;
            for (int i = 0; i < m_n_nodes; ++i)
                state.row(i) = init_state;
        }


        State operator()(const State &x,
                         const TArray2d &coupling,
                         const TArray1d &local_coupling) const override;

    };*/

}

#endif //TVB_CPP_MONTBRIO_H
