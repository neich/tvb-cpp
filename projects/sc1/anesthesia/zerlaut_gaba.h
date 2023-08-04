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

#ifndef TVB_CPP_ZERLAUT_GABA
#define TVB_CPP_ZERLAUT_GABA

#include <unsupported/Eigen/MatrixFunctions>

#include <tvb-cpp/definitions.h>
#include <tvb-cpp/simulator/models/zerlaut.h>

using namespace tvb;

class ZerlautGABA: public ZerlautAdaptationSecondOrder {
    TArray1d gaba_ratio;
//    = NArray(
//            label=r"gaba_ratio",
//    default=np.array([1.0]),
//    domain=Range(lo=0.0, hi=3.0, step=0.1),
//            doc="""Heterogeneous GABA density vector.""")

    TArray1d alpha_g;
//    = NArray(
//            label=r":math:`\alpha_g`",
//    default=np.array([1.0]),
//    domain=Range(lo=0.0, hi=3.0, step=0.1),
//            doc="""Linear component for heterogeneous tau_i generation.""")

    TArray1d beta_g;
//    = NArray(
//            label=r":math:`\beta_g`",
//    default=np.array([0.0]),
//    domain=Range(lo=-2.0, hi=2.0, step=0.1),
//            doc="""Constant component for heterogeneous tau_i generation.""")

    TArray1d gamma_g;
//    = NArray(
//            label=r":math:`\gamma_g`",
//    default=np.array([1.0]),
//    domain=Range(lo=0.0, hi=3.0, step=0.1),
//            doc="""Linear component for heterogeneous E_l generation.""")

    TArray1d delta_g;
//    = NArray(
//            label=r":math:`\delta_g`",
//    default=np.array([0.0]),
//    domain=Range(lo=-2.0, hi=2.0, step=0.1),
//            doc="""Constant component for heterogeneous E_l generation.""")

public:
    ZerlautGABA(int n_nodes) : ZerlautAdaptationSecondOrder(n_nodes) {
    }

    void configure() override {
        ZerlautAdaptationSecondOrder::configure();
        this->set_param("g_L", 10);
        this->set_param("E_L_e", -63.0);
        this->set_param("E_L_i", -65.0);
        this->set_param("C_m", 200);
        this->set_param("a_e", 0.0);
        this->set_param("b_e", 0.0);
        this->set_param("a_i", 0.0);
        this->set_param("b_i", 0.0);
        this->set_param("tau_w_e", 500.0);
        this->set_param("tau_w_i", 1.0);
        this->set_param("E_e", 0.0);
        this->set_param("E_i", -80.0);
        this->set_param("Q_e", 1.5);
        this->set_param("Q_i", 5.0);
        this->set_param("tau_e", 5.0);
        this->set_param("tau_i", 5.0);
        this->set_param("N_tot", 10000);
        this->set_param("p_connect_e", 0.05);
        this->set_param("p_connect_i", 0.05);
        this->set_param("g", 0.25);
        this->set_param("T", 20.0); // Changed from 40.0 to 20.0
        this->set_param("weight_noise", 1e-4);
        this->set_param("tau_OU", 5.0);

// TODO
//# I am changing this to put Fede's config, to see if it works better
//        'P_e':[-0.05017034,  0.00451531, -0.00794377, -0.00208418, -0.00054697,
//                0.00341614, -0.01156433,  0.00194753,  0.00274079, -0.01066769],
//        'P_i':[-0.05184978,  0.0061593 , -0.01403522,  0.00166511, -0.0020559 ,
//                0.00318432, -0.03112775,  0.00656668,  0.00171829, -0.04516385],
//        TArray1d P_e(10);
//        P_e << -0.0498,
//                0.00506,
//                -0.025,
//                0.0014,
//                -0.00041,
//                0.0105,
//                -0.036,
//                0.0074,
//                0.0012,
//                -0.0407;
//        TArray1d P_i(10);
//        P_i << -0.0514,
//                0.004,
//                -0.0083,
//                0.0002,
//                -0.0005,
//                0.0014,
//                -0.0146,
//                0.0045,
//                0.0028,
//                -0.0153;

        TArray1d P_e(10);
        P_e << -0.05017034,  0.00451531, -0.00794377, -0.00208418, -0.00054697,
                0.00341614, -0.01156433,  0.00194753,  0.00274079, -0.01066769;
        TArray1d P_i(10);
        P_i << -0.05184978,  0.0061593 , -0.01403522,  0.00166511, -0.0020559 ,
                0.00318432, -0.03112775,  0.00656668,  0.00171829, -0.04516385;

        this->set_param("P_e", P_e);
        this->set_param("P_i", P_i);
        this->set_param("K_ext_e", 400.0);
        this->set_param("K_ext_i", 0.0);

        this->set_param("external_input_ex_ex", 0.315*1e-3);
        this->set_param("external_input_ex_in", 0.0);
        this->set_param("external_input_in_ex", 0.315*1e-3);
        this->set_param("external_input_in_in", 0.0);

        gaba_ratio.resize(m_n_nodes);
        alpha_g.resize(m_n_nodes);
        beta_g.resize(m_n_nodes);
        delta_g.resize(m_n_nodes);
        gamma_g.resize(m_n_nodes);

        gaba_ratio.fill(1.0);
        alpha_g.fill(1.0);
        beta_g.fill(0.0);
        delta_g.fill(0.0);
        gamma_g.fill(1.0);
    }

    void set_param(const std::string& param, Float value) override {
        ADD_SETTER_FILL(gaba_ratio, alpha_g, beta_g, gamma_g, delta_g)
        ZerlautAdaptationSecondOrder::set_param(param, value);
    }

    void set_param(const std::string& param, const TArray1d& value) override {
        ADD_SETTER_VALUE(gaba_ratio, alpha_g, beta_g, gamma_g, delta_g)
        ZerlautAdaptationSecondOrder::set_param(param, value);
    }

    const TArray1d& get_param_value(const std::string& param) const override {
        ADD_GETTER(gaba_ratio, alpha_g, beta_g, gamma_g, delta_g)
        ZerlautAdaptationSecondOrder::get_param_value(param);
    }

    void init_dependant() override {
        this->tau_i = this->alpha_g * this->gaba_ratio + this->beta_g;
        this->E_L_e = this->gamma_g * this->gaba_ratio + this->delta_g;
    }

};


#endif //TVB_CPP_ZERLAUT_GABA
