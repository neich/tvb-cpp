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

#ifndef TVB_CPP_ZERLAUT
#define TVB_CPP_ZERLAUT

#include <unsupported/Eigen/MatrixFunctions>

#include <definitions.h>
#include <simulator/model.h>

namespace tvb {

    class ZerlautAdaptationFirstOrder : public Model {

//        r"""
//        **References**:
//        .. [ZD_2018]  Zerlaut, Y., Chemla, S., Chavane, F. et al. *Modeling mesoscopic cortical dynamics using a mean-field
//                model of conductance-based networks of adaptive
//        exponential integrate-and-fire neurons*,
//        J Comput Neurosci (2018) 44: 45. https://doi-org.lama.univ-amu.fr/10.1007/s10827-017-0668-2
//        .. [MV_2018]  Matteo di Volo, Alberto Romagnoni, Cristiano Capone, Alain Destexhe (2018)
//        *Mean-field model for the dynamics of conductance-based networks of excitatory and inhibitory spiking neurons
//                with adaptation*, bioRxiv, doi: https://doi.org/10.1101/352393
//
//        Used Eqns 4 from [MV_2018]_ in ``dfun``.
//
//        The default parameters are taken from table 1 of [ZD_2018]_, pag.47 and modify for the adaptation [MV_2018]
//        +---------------------------+------------+
//        |                 Table 1                |
//        +--------------+------------+------------+
//        |Parameter     |  Value     | Unit       |
//        +==============+============+============+
//        |             cellular property          |
//        +--------------+------------+------------+
//        | g_L          |   10.00    |   nS       |
//        +--------------+------------+------------+
//        | E_L_e        |  -60.00    |   mV       |
//        +--------------+------------+------------+
//        | E_L_i        |  -65.00    |   mV       |
//        +--------------+------------+------------+
//        | C_m          |   200.0    |   pF       |
//        +--------------+------------+------------+
//        | b_e          |   60.0     |   nS       |
//        +--------------+------------+------------+
//        | b_i          |   0.0      |   nS       |
//        +--------------+------------+------------+
//        | a_e          |   4.0      |   nS       |
//        +--------------+------------+------------+
//        | a_i          |   0.0      |   nS       |
//        +--------------+------------+------------+
//        | tau_w_e      |   500.0    |   ms       |
//        +--------------+------------+------------+
//        | tau_w_i      |   0.0      |   ms       |
//        +--------------+------------+------------+
//        | T            |   20.0      |   ms       |
//        +--------------+------------+------------+
//        |          synaptic properties           |
//        +--------------+------------+------------+
//        | E_e          |    0.0     | mV         |
//        +--------------+------------+------------+
//        | E_i          |   -80.0    | mV         |
//        +--------------+------------+------------+
//        | Q_e          |    1.0     | nS         |
//        +--------------+------------+------------+
//        | Q_i          |    5.0     | nS         |
//        +--------------+------------+------------+
//        | tau_e        |    5.0     | ms         |
//        +--------------+------------+------------+
//        | tau_i        |    5.0     | ms         |
//        +--------------+------------+------------+
//        |          numerical network             |
//        +--------------+------------+------------+
//        | N_tot        |  10000     |            |
//        +--------------+------------+------------+
//        | p_connect    |    5.0 %   |            |
//        +--------------+------------+------------+
//        | g            |   20.0 %   |            |
//        +--------------+------------+------------+
//        | K_e_ext      |   400      |            |
//        +--------------+------------+------------+
//        | K_i_ext      |   0        |            |
//        +--------------+------------+------------+
//        |external_input|    0.000   | Hz         |
//        +--------------+------------+------------+
//
//        The default coefficients of the transfer function are taken from table I of [MV_2018]_, pag.49
//        +-------------+-------------+-------------+-------------+-------------+-------------+-------------+-------------+-------------+-------------+
//        |      excitatory cell      |
//        +-------------+-------------+-------------+-------------+-------------+-------------+-------------+-------------+-------------+-------------+
//        |  -4.98e-02  |   5.06e-03  |  -2.5e-02   |   1.4e-03   |  -4.1e-04   |   1.05e-02  |  -3.6e-02   |   7.4e-03   |   1.2e-03   |  -4.07e-02  |
//        +-------------+-------------+-------------+-------------+-------------+-------------+-------------+-------------+-------------+-------------+
//        |      inhibitory cell      |
//        +-------------+-------------+-------------+-------------+-------------+-------------+-------------+-------------+-------------+-------------+
//        |  -5.14e-02  |   4.0e-03   |  -8.3e-03   |   2.0e-04   |  -5.0e-04   |   1.4e-03   |  -1.46e-02  |   4.5e-03   |   2.8e-03   |  -1.53e-02  |
//        +-------------+-------------+-------------+-------------+-------------+-------------+-------------+-------------+-------------+-------------+
//
//        The models (:math:`E`, :math:`I`) phase-plane, including a representation of
//        the vector field as well as its null-clines, using default parameters, can be
//        seen below:
//
//        .. automethod:: Zerlaut_adaptation_first_order.__init__
//
//                The general formulation for the Zerlaut adaptation first order model as a
//        dynamical unit at a node $k$ in a BNM with $l$ nodes reads:
//
//        .. math::
//        T\dot{E}_k &= F_e-E_k \\
//        T\dot{I}_k &= F_i-I_k \\
//        dot{W}_k &= W_k/tau_w_e-b*E_k \\
//        F_\lambda = Erfc(V^{eff}_{thre}-\mu_V/\sqrt(2)\sigma_V)
//
//        """

    public:
        TArray1d g_L ;
//        label=":math:`g_{L}`",
//        default=numpy.array([10.]),  # 10 nS by default, i.e. ~ 100MOhm input resitance at rest
//        domain=Range(lo=0.1, hi=100.0, step=0.1),  # 0.1nS would be a very small cell, 100nS a very big one
//        doc="""leak conductance [nS]""")

        TArray1d E_L_e;
//        label=":math:`E_{L}`",
//        default=numpy.array([-65.0]),
//        domain=Range(lo=-90.0, hi=-60.0, step=0.1),  # resting potential, usually between -85mV and -65mV
//        doc="""leak reversal potential for excitatory [mV]""")

        TArray1d E_L_i;
//        label=":math:`E_{L}`",
//        default=numpy.array([-65.0]),
//        domain=Range(lo=-90.0, hi=-60.0, step=0.1),  # resting potential, usually between -85mV and -65mV
//        doc="""leak reversal potential for inhibitory [mV]""")

        // N.B. Not independent of g_L, C_m should scale linearly with g_L
        TArray1d C_m;
//        = NArray(
//                label=":math:`C_{m}`",
//        default=numpy.array([200.0]),
//        domain=Range(lo=10.0, hi=500.0, step=10.0),  # 20pF very small cell, 400pF very
//                doc="""membrane capacitance [pF]""")

        TArray1d b_e;
//        label=":math:`b_e`",
//        default=numpy.array([60.0]),
//        domain=Range(lo=0.0, hi=150.0, step=1.0),
//                doc="""Excitatory adaptation current increment [pA]""")

        TArray1d a_e;
//        label=":math:`a_e`",
//        default=numpy.array([4.0]),
//        domain=Range(lo=0.0, hi=20.0, step=0.1),
//                doc="""Excitatory adaptation conductance [nS]""")

        TArray1d b_i;
//        label=":math:`b_i`",
//        default=numpy.array([0.0]),
//        domain=Range(lo=0.0, hi=100.0, step=0.1),
//                doc="""Inhibitory adaptation current increment [pA]""")

        TArray1d a_i;
//        label=":math:`a_i`",
//        default=numpy.array([0.0]),
//        domain=Range(lo=0.0, hi=20.0, step=0.1),
//                doc="""Inhibitory adaptation conductance [nS]""")

        TArray1d tau_w_e;
//        = NArray(
//                label=":math:`tau_w_e`",
//        default=numpy.array([500.0]),
//        domain=Range(lo=5.0, hi=1000.0, step=1.0),
//                doc="""Adaptation time constant [ms]""")

        TArray1d tau_w_i;
//        label=":math:`tau_wi`",
//        default=numpy.array([1.0]),
//        domain=Range(lo=1.0, hi=1000.0, step=1.0),
//                doc="""Adaptation time constant of inhibitory neurons [ms]""")

        TArray1d E_e;
//        = NArray(
//                label=r":math:`E_e`",
//        default=numpy.array([0.0]),
//        domain=Range(lo=-20., hi=20., step=0.01),
//                doc="""excitatory reversal potential [mV]""")

        TArray1d E_i;
//        = NArray(
//                label=":math:`E_i`",
//        default=numpy.array([-80.0]),
//        domain=Range(lo=-100.0, hi=-60.0, step=1.0),
//                doc="""inhibitory reversal potential [mV]""")

        TArray1d Q_e;
//        = NArray(
//                label=r":math:`Q_e`",
//        default=numpy.array([1.5]),
//        domain=Range(lo=0.0, hi=5.0, step=0.1),
//                doc="""excitatory quantal conductance [nS]""")

        TArray1d Q_i;
//        = NArray(
//                label=r":math:`Q_i`",
//        default=numpy.array([5.0]),
//        domain=Range(lo=0.0, hi=10.0, step=0.1),
//                doc="""inhibitory quantal conductance [nS]""")

        TArray1d tau_e;
//        = NArray(
//                label=":math:`\tau_e`",
//        default=numpy.array([5.0]),
//        domain=Range(lo=1.0, hi=10.0, step=1.0),
//                doc="""excitatory decay [ms]""")

        TArray1d tau_i;
//        = NArray(
//                label=":math:`\tau_i`",
//        default=numpy.array([5.0]),
//        domain=Range(lo=0.5, hi=10.0, step=0.01),
//                doc="""inhibitory decay [ms]""")

        TArray1d N_tot;
//        = NArray(
//                dtype=numpy.int,
//        label=":math:`N_{tot}`",
//        default=numpy.array([10000]),
//        domain=Range(lo=1000, hi=50000, step=1000),
//                doc="""cell number""")

        TArray1d p_connect;
//        = NArray(
//                label=":math:`\epsilon`",
//        default=numpy.array([0.05]),
//        domain=Range(lo=0.001, hi=0.2, step=0.001),  # valid only for relatively sparse connectivities
//                doc="""connectivity probability""")

        TArray1d g;
//        = NArray(
//                label=":math:`g`",
//        default=numpy.array([0.2]),
//        domain=Range(lo=0.01, hi=0.4, step=0.01),  # inhibitory cell number never overcomes excitatory ones
//                doc="""fraction of inhibitory cells""")

        TArray1d K_ext_e;
//        label=":math:`K_ext_e`",
//        default=numpy.array([400]),
//        domain=Range(lo=0, hi=10000, step=1),  # inhibitory cell number never overcomes excitatory ones
//                doc="""Number of excitatory connexions from external population""")

        TArray1d K_ext_i;
//        label=":math:`K_ext_i`",
//        default=numpy.array([0]),
//        domain=Range(lo=0, hi=10000, step=1),  # inhibitory cell number never overcomes excitatory ones
//                doc="""Number of inhibitory connexions from external population""")

        TArray1d T;
//        = NArray(
//                label=":math:`T`",
//        default=numpy.array([20.0]),
//        domain=Range(lo=1., hi=20.0, step=0.1),
//                doc="""time scale of describing network activity""")

        TArray1d P_e;
//        = NArray(
//                label=":math:`P_e`",  # TODO need to check the size of the array when it's used
//        default=numpy.array([-4.98e-02, 5.06e-03, -2.5e-02, 1.4e-03,
//        -4.1e-04, 1.05e-02, -3.6e-02, 7.4e-03,
//        1.2e-03, -4.07e-02]),
//        doc="""Polynome of excitatory phenomenological threshold (order 9)""")

        TArray1d P_i;
//        = NArray(
//                label=":math:`P_i`",  # TODO need to check the size of the array when it's used
//        default=numpy.array([-5.14e-02, 4.0e-03, -8.3e-03, 2.0e-04,
//        -5.0e-04, 1.4e-03, -1.46e-02, 4.5e-03,
//        2.8e-03, -1.53e-02]),
//        doc="""Polynome of inhibitory phenomenological threshold (order 9)""")

        TArray1d external_input_ex_ex;
//        label=":math:`\nu_e^{drive}`",
//        default=numpy.array([0.000]),
//        domain=Range(lo=0.00, hi=0.1, step=0.001),
//                doc="""external drive""")

        TArray1d external_input_ex_in;
//        label=":math:`\nu_e^{drive}`",
//        default=numpy.array([0.000]),
//        domain=Range(lo=0.00, hi=0.1, step=0.001),
//                doc="""external drive""")

        TArray1d external_input_in_ex;
//        label=":math:`\nu_e^{drive}`",
//        default=numpy.array([0.000]),
//        domain=Range(lo=0.00, hi=0.1, step=0.001),
//                doc="""external drive""")

        TArray1d external_input_in_in;
//        label=":math:`\nu_e^{drive}`",
//        default=numpy.array([0.000]),
//        domain=Range(lo=0.00, hi=0.1, step=0.001),
//                doc="""external drive""")



    public:
        explicit ZerlautAdaptationFirstOrder(int n_nodes) : Model(n_nodes, 4) {
            this->configure(n_nodes);
            m_n_vars = 4;
            m_cvars = { 0 };
            m_state_vars = { "E", "I", "W_e", "W_i" };
        }

        void set_param(const std::string& param, Float value) {
            ADD_SETTER(g_L, E_L_e, E_L_i, C_m, b_e, a_e, b_i, a_i, tau_w_e, tau_w_i,\
                       E_e, E_i, Q_e, Q_i, tau_e, tau_i, N_tot, p_connect, g, K_ext_e,\
                       K_ext_i, T, external_input_ex_ex, external_input_ex_in,\
                       external_input_in_ex, external_input_in_in)
        }

        void configure(int n_nodes) {
            g_L.resize(n_nodes);
            E_L_e.resize(n_nodes);
            E_L_i.resize(n_nodes);
            C_m.resize(n_nodes);
            b_e.resize(n_nodes);
            a_e.resize(n_nodes);
            b_i.resize(n_nodes);
            a_i.resize(n_nodes);
            tau_w_e.resize(n_nodes);
            tau_w_i.resize(n_nodes);
            E_e.resize(n_nodes);
            E_i.resize(n_nodes);
            Q_e.resize(n_nodes);
            Q_i.resize(n_nodes);
            tau_e.resize(n_nodes);
            tau_i.resize(n_nodes);
            N_tot.resize(n_nodes);
            p_connect.resize(n_nodes);
            g.resize(n_nodes);
            K_ext_e.resize(n_nodes);
            K_ext_i.resize(n_nodes);
            T.resize(n_nodes);
            P_e.resize(10);
            P_i.resize(10);
            external_input_ex_ex.resize(n_nodes);
            external_input_ex_in.resize(n_nodes);
            external_input_in_ex.resize(n_nodes);
            external_input_in_in.resize(n_nodes);


            g_L.fill(10.0);
            E_L_e.fill(-65.0);
            E_L_i.fill(-65.0);
            C_m.fill(200.0);
            b_e.fill(60.0);
            a_e.fill(4.0);
            b_i.fill(0.0);
            a_i.fill(0.0);
            tau_w_e.fill(500.0);
            tau_w_i.fill(1.0);
            E_e.fill(0.0);
            E_i.fill(-80.0);
            Q_e.fill(1.5);
            Q_i.fill(5.0);
            tau_e.fill(5.0);
            tau_i.fill(5.0);
            N_tot.fill(10000);
            p_connect.fill(0.05);
            g.fill(0.2);
            K_ext_e.fill(400);
            K_ext_i.fill(0);
            T.fill(20.0);
            P_e << -0.04983106, 0.005063550882777035, -0.023470121807314552,
                    0.0022951513725067503,
                    -0.0004105302652029825, 0.010547051343547399, -0.03659252821136933,
                    0.007437487505797858, 0.001265064721846073, -0.04072161294490446;
            P_i << -0.05149122024209484, 0.004003689190271077, -0.008352013668528155,
                    0.0002414237992765705,
                    -0.0005070645080016026, 0.0014345394104282397, -0.014686689498949967,
                    0.004502706285435741,
                    0.0028472190352532454, -0.015357804594594548;
            external_input_ex_ex.fill(0.0);
            external_input_ex_in.fill(0.0);
            external_input_in_ex.fill(0.0);
            external_input_in_in.fill(0.0);
        }

        [[nodiscard]] State initial() const override {
            State result(m_n_nodes, m_n_vars);
            this->initial(result);
            return result;
        }

        void initial(State& state) const override {
            state.resize(m_n_nodes, m_n_vars);
            TArray1d init_state(m_n_vars);
            init_state << 0.1, 0.1, 100.0, 0.0;
            for (int i = 0; i < m_n_nodes; ++i)
                state.row(i) = init_state;
        }



        State operator()(const State &x,
                const TArray2d &coupling,
                const TArray1d &local_coupling) const override;

        inline TArray1d TF_excitatory(const TArray1d& fe, const TArray1d& fi, const TArray1d& fe_ext, const TArray1d& fi_ext, const TArray1d& W) const;

        inline TArray1d TF_inhibitory(const TArray1d& fe, const TArray1d& fi, const TArray1d& fe_ext, const TArray1d& fi_ext, const TArray1d& W) const;

        [[nodiscard]] TArray1d TF(const TArray1d& fe, const TArray1d& fi, const TArray1d& fe_ext, const TArray1d& fi_ext, const TArray1d& W, const TArray1d &P, const TArray1d &E_L) const;
    };

    class ZerlautAdptationSecondOrder: public ZerlautAdaptationFirstOrder {
    public:
        ZerlautAdptationSecondOrder(int n_nodes) : ZerlautAdaptationFirstOrder(n_nodes) {
            this->configure(n_nodes);
            m_n_vars = 7;
            m_cvars = { 0 };
            m_state_vars = { "E", "I", "C_ee", "C_ei", "C_ii", "W_e", "W_i" };
        }

        void initial(State& state) const override {
            state.resize(m_n_nodes, m_n_vars);
            TArray1d init_state(m_n_vars);
            init_state << 0.05, 0.05, 0.0, 0.0, 0.0, 100.0, 0.0;
            for (int i = 0; i < m_n_nodes; ++i)
                state.row(i) = init_state;
        }

        State operator()(const State &x,
                         const TArray2d &coupling,
                         const TArray1d &local_coupling) const override;

    protected:
        // Derivatives taken numerically : use a central difference formula with spacing `dx`
        inline TArray1d _diff_fe_E(const TArray1d& fe, const TArray1d& fi, const TArray1d& fe_ext, const TArray1d& fi_ext, const TArray1d& W, double df=1e-7) const {
            return (this->TF_excitatory(fe + df, fi, fe_ext, fi_ext, W) - this->TF_excitatory(fe - df, fi, fe_ext, fi_ext, W)) / (2 * df * 1e3);
        }
        inline TArray1d _diff_fe_I(const TArray1d& fe, const TArray1d& fi, const TArray1d& fe_ext, const TArray1d& fi_ext, const TArray1d& W, double df=1e-7) const {
            return (this->TF_inhibitory(fe + df, fi, fe_ext, fi_ext, W) - this->TF_inhibitory(fe - df, fi, fe_ext, fi_ext, W)) / (2 * df * 1e3);
        }

        inline TArray1d _diff_fi_E(const TArray1d& fe, const TArray1d& fi, const TArray1d& fe_ext, const TArray1d& fi_ext, const TArray1d& W, double df=1e-7) const {
            return (this->TF_excitatory(fe, fi + df, fe_ext, fi_ext, W) - this->TF_excitatory(fe, fi - df, fe_ext, fi_ext, W)) / (2 * df * 1e3);
        }
        inline TArray1d _diff_fi_I(const TArray1d& fe, const TArray1d& fi, const TArray1d& fe_ext, const TArray1d& fi_ext, const TArray1d& W, double df=1e-7) const {
            return (this->TF_inhibitory(fe, fi + df, fe_ext, fi_ext, W) - this->TF_inhibitory(fe, fi - df, fe_ext, fi_ext, W)) / (2 * df * 1e3);
        }

        inline TArray1d _diff2_fe_fe_e(const TArray1d& fe, const TArray1d& fi, const TArray1d& fe_ext, const TArray1d& fi_ext, const TArray1d& W, const TArray1d& _TF, double df=1e-7) const {
            return (this->TF_excitatory(fe + df, fi, fe_ext, fi_ext, W) - 2 * _TF + this->TF_excitatory(fe - df, fi, fe_ext, fi_ext, W)) /
                   ((df * 1e3) * (df * 1e3));
        }

        inline TArray1d _diff2_fe_fe_i(const TArray1d& fe, const TArray1d& fi, const TArray1d& fe_ext, const TArray1d& fi_ext, const TArray1d& W, const TArray1d& _TF, double df=1e-7) const {
            return (this->TF_inhibitory(fe + df, fi, fe_ext, fi_ext, W) - 2 * _TF + this->TF_inhibitory(fe - df, fi, fe_ext, fi_ext, W)) /
                   ((df * 1e3) * (df * 1e3));
        }

        inline TArray1d _diff2_fi_fe_E(const TArray1d& fe, const TArray1d& fi, const TArray1d& fe_ext, const TArray1d& fi_ext, const TArray1d& W, double df=1e-7) const {
            return (_diff_fi_E(fe + df, fi, fe_ext, fi_ext, W) - _diff_fi_E(fe - df, fi, fe_ext, fi_ext, W)) / (2 * df * 1e3);
        }
        inline TArray1d _diff2_fi_fe_I(const TArray1d& fe, const TArray1d& fi, const TArray1d& fe_ext, const TArray1d& fi_ext, const TArray1d& W, double df=1e-7) const {
            return (_diff_fi_I(fe + df, fi, fe_ext, fi_ext, W) - _diff_fi_I(fe - df, fi, fe_ext, fi_ext, W)) / (2 * df * 1e3);
        }

        inline TArray1d _diff2_fe_fi_E(const TArray1d& fe, const TArray1d& fi, const TArray1d& fe_ext, const TArray1d& fi_ext, const TArray1d& W, double df=1e-7) const {
            return (_diff_fe_E(fe, fi+df, fe_ext, fi_ext, W)-_diff_fe_E(fe, fi-df, fe_ext, fi_ext, W))/(2*df*1e3);
        }
        inline TArray1d _diff2_fe_fi_I(const TArray1d& fe, const TArray1d& fi, const TArray1d& fe_ext, const TArray1d& fi_ext, const TArray1d& W, double df=1e-7) const {
            return (_diff_fe_I(fe, fi+df, fe_ext, fi_ext, W)-_diff_fe_I(fe, fi-df, fe_ext, fi_ext, W))/(2*df*1e3);
        }

        inline TArray1d _diff2_fi_fi_e(const TArray1d& fe, const TArray1d& fi, const TArray1d& fe_ext, const TArray1d& fi_ext, const TArray1d& W, const TArray1d& _TF, double df=1e-7) const {
            return (this->TF_excitatory(fe, fi + df, fe_ext, fi_ext, W) - 2 * _TF + this->TF_excitatory(fe, fi - df, fe_ext, fi_ext, W)) /
                   ((df * 1e3) * (df * 1e3));
        }

        inline TArray1d _diff2_fi_fi_i(const TArray1d& fe, const TArray1d& fi, const TArray1d& fe_ext, const TArray1d& fi_ext, const TArray1d& W, const TArray1d& _TF, double df=1e-7) const {
            return (this->TF_inhibitory(fe, fi + df, fe_ext, fi_ext, W) - 2 * _TF + this->TF_inhibitory(fe, fi - df, fe_ext, fi_ext, W)) /
                   ((df * 1e3) * (df * 1e3));
        }
    };
}

#endif //TVB_CPP_ZERLAUT
