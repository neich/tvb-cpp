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

#include <cmath>

#include <simulator/models/zerlaut.h>

using namespace tvb;

std::tuple<TArray1d, TArray1d, TArray1d> get_fluct_regime_vars(const TArray1d &Fe,
                                                               const TArray1d &Fi,
                                                               const TArray1d &Fe_ext,
                                                               const TArray1d &Fi_ext,
                                                               const TArray1d &W,
                                                               const TArray1d &Q_e,
                                                               const TArray1d &tau_e,
                                                               const TArray1d &E_e,
                                                               const TArray1d &Q_i,
                                                               const TArray1d &tau_i,
                                                               const TArray1d &E_i,
                                                               const TArray1d &g_L,
                                                               const TArray1d &C_m,
                                                               const TArray1d &E_L,
                                                               const TArray1d &N_tot,
                                                               const TArray1d &p_connect,
                                                               const TArray1d &g,
                                                               const TArray1d &K_ext_e,
                                                               const TArray1d &K_ext_i);

TArray1d threshold_func(const TArray1d &muV,
                        const TArray1d &sigmaV,
                        const TArray1d &TvN,
                        double P0, double P1, double P2, double P3, double P4, double P5,
                        double P6, double P7, double P8, double P9);

TArray1d estimate_firing_rate(const TArray1d &muV, const TArray1d &sigmaV,
                              const TArray1d &Tv, const TArray1d &Vthre);


State ZerlautAdaptationFirstOrder::operator()(const State &x,
                                              const TArray2d &coupling,
                                              const TArray1d &local_coupling) {

    State derivative(m_n_nodes, m_n_vars);

    const TArray1d &E = x.col(0);
    const TArray1d &I = x.col(1);
    const TArray1d &W_e = x.col(2);
    const TArray1d &W_i = x.col(3);

    const TArray1d &c_0 = coupling.col(0);

    TArray1d lc_E = local_coupling * E;
    TArray1d lc_I = local_coupling * I;

    TArray1d Fe_ext = c_0 + lc_E;
    TArray1d Fi_ext = lc_I;

    // Excitatory firing rate derivation
    TArray1d tmp = (
            this->TF_excitatory(E, I, Fe_ext + this->external_input_ex_ex, Fi_ext + this->external_input_ex_in, W_e) -
            E);
    derivative.col(0) = tmp / this->T;
    // Inhibitory firing rate derivation
    tmp = (this->TF_inhibitory(E, I, Fe_ext + this->external_input_in_ex, Fi_ext + this->external_input_in_in, W_i) -
           I);
    derivative.col(1) = tmp / this->T;
    // Adaptation excitatory
    auto [mu_V, sigma_V, T_V] = get_fluct_regime_vars(
            E, I, Fe_ext + this->external_input_ex_ex,
            Fi_ext + this->external_input_ex_in,
            W_e, this->Q_e, this->tau_e, this->E_e,
            this->Q_i, this->tau_i, this->E_i,
            this->g_L, this->C_m, this->E_L_e, this->N_tot,
            this->p_connect, this->g, this->K_ext_e, this->K_ext_i);
    derivative.col(2) = -W_e / this->tau_w_e + this->b_e * E + this->a_e * (mu_V - this->E_L_e) / this->tau_w_e;
    // Adaptation inhibitory
    auto [mu_V_i, sigma_V_i, T_V_i] = get_fluct_regime_vars(
            E, I, Fe_ext + this->external_input_in_ex,
            Fi_ext + this->external_input_in_in,
            W_i, this->Q_e, this->tau_e, this->E_e,
            this->Q_i, this->tau_i, this->E_i,
            this->g_L, this->C_m, this->E_L_i, this->N_tot,
            this->p_connect, this->g, this->K_ext_e, this->K_ext_i);
    derivative.col(3) = -W_i / this->tau_w_i + this->b_i * I + this->a_i * (mu_V_i - this->E_L_i) / this->tau_w_i;


    return derivative;

}

inline
TArray1d ZerlautAdaptationFirstOrder::TF_excitatory(const TArray1d &fe, const TArray1d &fi, const TArray1d &fe_ext,
                                                    const TArray1d &fi_ext, const TArray1d &W) const {
    return this->TF(fe, fi, fe_ext, fi_ext, W, this->P_e, this->E_L_e);
}

inline
TArray1d ZerlautAdaptationFirstOrder::TF_inhibitory(const TArray1d &fe, const TArray1d &fi, const TArray1d &fe_ext,
                                                    const TArray1d &fi_ext, const TArray1d &W) const {
    return this->TF(fe, fi, fe_ext, fi_ext, W, this->P_i, this->E_L_i);
}

inline
TArray1d
ZerlautAdaptationFirstOrder::TF(const TArray1d &fe, const TArray1d &fi, const TArray1d &fe_ext, const TArray1d &fi_ext,
                                const TArray1d &W,
                                const TArray1d &P, const TArray1d &E_L) const {

    //transfer function for inhibitory population
    //Inspired from the next repository :
    //https://github.com/yzerlaut/notebook_papers/tree/master/modeling_mesoscopic_dynamics
    //:param fe: firing rate of excitatory population
    //:param fi: firing rate of inhibitory population
    //:param W: level of adaptation
    //:param P: Polynome of neurons phenomenological threshold (order 9)
    //:param E_L: leak reversal potential
    //:return: result of transfer function

    auto [mu_V, sigma_V, T_V] =
            get_fluct_regime_vars(fe, fi, fe_ext, fi_ext, W, this->Q_e, this->tau_e, this->E_e,
                                  this->Q_i, this->tau_i, this->E_i,
                                  this->g_L, this->C_m, E_L, this->N_tot,
                                  this->p_connect, this->g, this->K_ext_e, this->K_ext_i);
    TArray1d V_thre = threshold_func(mu_V, sigma_V, (T_V * this->g_L) / this->C_m,
                                     P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7], P[8], P[9]);
    V_thre *= 1e3;  // the threshold need to be in mv and not in Volt
    TArray1d f_out = estimate_firing_rate(mu_V, sigma_V, T_V, V_thre);
    return f_out;
}

std::tuple<TArray1d, TArray1d, TArray1d> get_fluct_regime_vars(const TArray1d &Fe,
                                                               const TArray1d &Fi,
                                                               const TArray1d &Fe_ext,
                                                               const TArray1d &Fi_ext,
                                                               const TArray1d &W,
                                                               const TArray1d &Q_e,
                                                               const TArray1d &tau_e,
                                                               const TArray1d &E_e,
                                                               const TArray1d &Q_i,
                                                               const TArray1d &tau_i,
                                                               const TArray1d &E_i,
                                                               const TArray1d &g_L,
                                                               const TArray1d &C_m,
                                                               const TArray1d &E_L,
                                                               const TArray1d &N_tot,
                                                               const TArray1d &p_connect,
                                                               const TArray1d &g,
                                                               const TArray1d &K_ext_e,
                                                               const TArray1d &K_ext_i) {
    //"""
    //Compute the mean characteristic of neurons.
    //Inspired from the next repository :
    //https://github.com/yzerlaut/notebook_papers/tree/master/modeling_mesoscopic_dynamics
    //:param Fe: firing rate of excitatory population
    //:param Fi: firing rate of inhibitory population
    //:param W: level of adaptation
    //:param Q_e: excitatory quantal conductance
    //:param tau_e: excitatory decay
    //:param E_e: excitatory reversal potential
    //:param Q_i: inhibitory quantal conductance
    //:param tau_i: inhibitory decay
    //:param E_i: inhibitory reversal potential
    //:param E_L: leakage reversal voltage of neurons
    //:param g_L: leak conductance
    //:param C_m: membrane capacitance
    //:param E_L: leak reversal potential
    //:param N_tot: cell number
    //:param p_connect: connectivity probability
    //:param g: fraction of inhibitory cells
    //:return: mean and variance of membrane voltage of neurons and autocorrelation time constant
    //"""
    // firing rate
    // 1e-6 represent spontaneous release of synaptic neurotransmitter or some intrinsic currents of neurons
    TArray1d fe = (Fe + 1e-6) * (1. - g) * p_connect * N_tot + Fe_ext * K_ext_e;
    TArray1d fi = (Fi + 1e-6) * g * p_connect * N_tot + Fi_ext * K_ext_i;

    // conductance fluctuation and effective membrane time constant
    TArray1d mu_Ge = Q_e * tau_e * fe;
    TArray1d mu_Gi = Q_i * tau_i * fi;  // Eqns 5 from [MV_2018]
    TArray1d mu_G = g_L + mu_Ge + mu_Gi;  // Eqns 6 from [MV_2018]
    TArray1d T_m = C_m / mu_G; // Eqns 6 from [MV_2018]

    // membrane potential
    TArray1d mu_V = (mu_Ge * E_e + mu_Gi * E_i + g_L * E_L - W) / mu_G;  // Eqns 7 from [MV_2018]
    // post-synaptic membrane potential event s around muV
    TArray1d U_e = Q_e / mu_G * (E_e - mu_V);
    TArray1d U_i = Q_i / mu_G * (E_i - mu_V);
    // Standard deviation of the fluctuations
    // Eqns 8 from [MV_2018]
    TArray1d sigma_V =
            fe * (U_e * tau_e).pow(2.0) / (2. * (tau_e + T_m)) + fi * (U_i * tau_i).pow(2.0) / (2. * (tau_i + T_m));
    sigma_V = sigma_V.sqrt();
    // Autocorrelation-time of the fluctuations Eqns 9 from [MV_2018]
    TArray1d T_V_numerator = (fe * (U_e * tau_e).pow(2.0) + fi * (U_i * tau_i).pow(2.0));
    TArray1d T_V_denominator = (fe * (U_e * tau_e).pow(2.0) / (tau_e + T_m) +
                                fi * (U_i * tau_i).pow(2.0) / (tau_i + T_m));
    TArray1d T_V = T_V_numerator.binaryExpr(T_V_denominator, [](Float n, Float d) { return d != 0.0 ? n / d : n; });
    //TArray1d T_V = T_V_numerator.binaryExpr(T_V_denominator, [](double n, double d) { return d != 0.0 ? n/d : n; });
    return {mu_V, sigma_V, T_V};
}

TArray1d threshold_func(const TArray1d &muV,
                        const TArray1d &sigmaV,
                        const TArray1d &TvN,
                        double P0, double P1, double P2, double P3, double P4, double P5,
                        double P6, double P7, double P8, double P9) {
    //The threshold function of the neurons
    //:param muV: mean of membrane voltage
    //:param sigmaV: variance of membrane voltage
    //:param TvN: autocorrelation time constant
    //:param P: Fitted coefficients of the transfer functions
    //:return: threshold of neurons

    // Normalization factors page 48 after the equation 4 from [ZD_2018]
    double muV0 = -60.0;
    double DmuV0 = 10.0;
    double sV0 = 4.0;
    double DsV0 = 6.0;
    double TvN0 = 0.5;
    double DTvN0 = 1.0;
    TArray1d V = (muV - muV0) / DmuV0;
    TArray1d S = (sigmaV - sV0) / DsV0;
    TArray1d T = (TvN - TvN0) / DTvN0;
    // Eqns 11 from [MV_2018]
    return P0 + P1 * V + P2 * S + P3 * T + P4 * V.pow(2.0) + P5 * S.pow(2.0) + P6 * T.pow(2.0) + P7 * V * S +
           P8 * V * T +
           P9 * S * T;
}

TArray1d estimate_firing_rate(const TArray1d &muV, const TArray1d &sigmaV,
                              const TArray1d &Tv, const TArray1d &Vthre) {
    // Eqns 10 from [MV_2018]
    TArray1d e = ((Vthre - muV) / (M_SQRT2 * sigmaV)).unaryExpr([](Float x) { return std::erfc(x); });
    // TArray1d e = ((Vthre - muV) / (M_SQRT2 * sigmaV)).unaryExpr([](double x){return std::erfc(x);});
    return e / (Tv * 2.0);
}

State ZerlautAdptationSecondOrder::operator()(const State &x, const TArray2d &coupling,
                                              const TArray1d &local_coupling) {
//    .. math::
//    \forall \mu,\lambda,\eta \in \{e,i\}^3\,
//    \left\{
//        \begin{split}
//        T \, \frac{\partial \nu_\mu}{\partial t} = & (\mathcal{F}_\mu - \nu_\mu )
//        + \frac{1}{2} \, c_{\lambda \eta} \,
//        \frac{\partial^2 \mathcal{F}_\mu}{\partial \nu_\lambda \partial \nu_\eta} \\
//            T \, \frac{\partial c_{\lambda \eta} }{\partial t}  =  & A_{\lambda \eta} +
//        (\mathcal{F}_\lambda - \nu_\lambda ) \, (\mathcal{F}_\eta - \nu_\eta ) + \\
//            & c_{\lambda \mu} \frac{\partial \mathcal{F}_\mu}{\partial \nu_\lambda} +
//                c_{\mu \eta} \frac{\partial \mathcal{F}_\mu}{\partial \nu_\eta}
//        - 2  c_{\lambda \eta}
//        \end{split}
//        \right.
//                dot{W}_k &= W_k/tau_w-b*E_k  \\
//
//        with:
//        A_{\lambda \eta} =
//        \left\{
//            \begin{split}
//            \frac{\mathcal{F}_\lambda \, (1/T - \mathcal{F}_\lambda)}{N_\lambda}
//            \qquad & \textrm{if  } \lambda=\eta \\
//            0 \qquad & \textrm{otherwise}
//            \end{split}
//            \right.


    TArray1d N_e = this->N_tot * (1.0 - this->g);
    TArray1d N_i = this->N_tot * this->g;

    State derivative(m_n_nodes, m_n_vars);

    const TArray1d &E = x.col(0);
    const TArray1d &I = x.col(1);
    const TArray1d &C_ee = x.col(2);
    const TArray1d &C_ei = x.col(3);
    const TArray1d &C_ii = x.col(4);
    const TArray1d &W_e = x.col(5);
    const TArray1d &W_i = x.col(6);
//
//// long-range coupling
    const TArray1d &c_0 = coupling.col(0);
//
//// short-range (local) coupling
    TArray1d lc_E = local_coupling * E;
    TArray1d lc_I = local_coupling * I;
//
    TArray1d E_input_excitatory = c_0 + lc_E + this->external_input_ex_ex;
    TArray1d E_input_inhibitory = c_0 + lc_E + this->external_input_in_ex;
    TArray1d I_input_excitatory = lc_I + this->external_input_ex_in;
    TArray1d I_input_inhibitory = lc_I + this->external_input_in_in;
//
//// Transfer function of excitatory and inhibitory neurons
    TArray1d TF_e = this->TF_excitatory(E, I, E_input_excitatory, I_input_excitatory, W_e);
    TArray1d TF_i = this->TF_inhibitory(E, I, E_input_inhibitory, I_input_inhibitory, W_i);
//
//
//// Precompute some result
    TArray1d _diff_fe_TF_e = this->_diff_fe_E(E, I, E_input_excitatory, I_input_excitatory, W_e);
    TArray1d _diff_fe_TF_i = this->_diff_fe_I(E, I, E_input_inhibitory, I_input_inhibitory, W_i);
    TArray1d _diff_fi_TF_e = this->_diff_fi_E(E, I, E_input_excitatory, I_input_excitatory, W_e);
    TArray1d _diff_fi_TF_i = this->_diff_fi_I(E, I, E_input_inhibitory, I_input_inhibitory, W_i);
//
//// equation is inspired from github of Zerlaut :
//// https://github.com/yzerlaut/notebook_papers/blob/master/modeling_mesoscopic_dynamics/mean_field/master_equation.py
//// Excitatory firing rate derivation
    derivative.col(0) = (TF_e - E
                         + .5 * C_ee * this->_diff2_fe_fe_e(E, I, E_input_excitatory, I_input_excitatory, W_e, TF_e)
                         + .5 * C_ei * this->_diff2_fe_fi_E(E, I, E_input_excitatory, I_input_excitatory, W_e)
                         + .5 * C_ei * this->_diff2_fi_fe_E(E, I, E_input_excitatory, I_input_excitatory, W_e)
                         + .5 * C_ii * this->_diff2_fi_fi_e(E, I, E_input_excitatory, I_input_excitatory, W_e, TF_e)
                        ) / this->T;
//// Inhibitory firing rate derivation
    derivative.col(1) = (TF_i - I
                         + .5 * C_ee * this->_diff2_fe_fe_i(E, I, E_input_inhibitory, I_input_inhibitory, W_i, TF_i)
                         + .5 * C_ei * this->_diff2_fe_fi_I(E, I, E_input_inhibitory, I_input_inhibitory, W_i)
                         + .5 * C_ei * this->_diff2_fi_fe_I(E, I, E_input_inhibitory, I_input_inhibitory, W_i)
                         + .5 * C_ii * this->_diff2_fi_fi_i(E, I, E_input_inhibitory, I_input_inhibitory, W_i, TF_i)
                        ) / this->T;
//// Covariance excitatory-excitatory derivation
    derivative.col(2) = (TF_e * (1. / this->T - TF_e) / N_e
                         + (TF_e - E).pow(2.0)
                         + 2. * C_ee * _diff_fe_TF_e
                         + 2. * C_ei * _diff_fi_TF_i
                         - 2. * C_ee
                        ) / this->T;
//// Covariance excitatory-inhibitory or inhibitory-excitatory derivation
    derivative.col(3) = ((TF_e - E) * (TF_i - I)
                         + C_ee * _diff_fe_TF_e
                         + C_ei * _diff_fe_TF_i
                         + C_ei * _diff_fi_TF_e
                         + C_ii * _diff_fi_TF_i
                         - 2. * C_ei
                        ) / this->T;
//// Covariance inhibitory-inhibitory derivation
    derivative.col(4) = (TF_i * (1. / this->T - TF_i) / N_i
                         + (TF_i - I).pow(2.0)
                         + 2. * C_ii * _diff_fi_TF_i
                         + 2. * C_ei * _diff_fe_TF_e
                         - 2. * C_ii
                        ) / this->T;
//// Adaptation excitatory
    auto [mu_V, sigma_V, T_V] = get_fluct_regime_vars(
            E, I,
            E_input_excitatory,
            E_input_inhibitory,
            W_e, this->Q_e, this->tau_e, this->E_e,
            this->Q_i, this->tau_i, this->E_i,
            this->g_L, this->C_m, this->E_L_e, this->N_tot,
            this->p_connect, this->g, this->K_ext_e, this->K_ext_i);
    derivative.col(5) = -W_e / this->tau_w_e + this->b_e * E + this->a_e * (mu_V - this->E_L_e) / this->tau_w_e;

    //// Adaptation inhibitory
    auto [mu_V_i, sigma_V_i, T_V_i] = get_fluct_regime_vars(
            E, I,
            I_input_excitatory,
            I_input_inhibitory,
            W_i, this->Q_e, this->tau_e, this->E_e,
            this->Q_i, this->tau_i, this->E_i,
            this->g_L, this->C_m, this->E_L_i, this->N_tot,
            this->p_connect, this->g, this->K_ext_e, this->K_ext_i);
    derivative.col(6) = -W_i / this->tau_w_i + this->b_i * I + this->a_i * (mu_V_i - this->E_L_i) / this->tau_w_i;

    return derivative;
}
