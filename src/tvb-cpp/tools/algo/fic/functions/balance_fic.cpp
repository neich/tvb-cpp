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

#include "balance_fic.h"

#include <random>


using namespace tvb;

std::tuple<int, double, TArray1d> optimize_fic(SimConfig &sim_config, int voi, float value_base) {

    int N = sim_config.connectivity()->weights().rows();

    auto *model = sim_config.model();

    TArray1d delta(N);
    delta.fill(1.0);
    TArray1d distance(N);
    distance.fill(10.0);
    State initial_state;
    double prev_largest_distance = 0;
    // double sf = 1.0;
    double speed = 1.0;

    TArray1d J_i = model->get_param_value("J_i");
    TArray1d prev_J_i = J_i;
    double best_distance = 1e6;
    TArray1d best_Ji;
    OptResult result;

    double slow_factor = 1.0;

    unsigned step = 0;
    for (; step < 500; ++step) {
        //for (unsigned step = 0; step < 500; ++step) {
        model->set_param("J_i", J_i);
        auto [converged, sim_result] = tvb::simulate(sim_config, 1.0, voi);

        unsigned num_above_error = 0;
        unsigned n_states = sim_result->getRecords().size();
        double largest_distance = 0.0;
        unsigned largest_distance_i = 0;
        TArray1d distance(N);
        TArray1d min_distance(N);

        TArray1d total_ie = TArray1d::Zero(N);
        int skip_n_states = int(0.8 * n_states);
        for (unsigned j = skip_n_states; j < n_states; ++j)
            total_ie += sim_result->getRecords()[j].record.col(0);
        TArray1d ie = total_ie / (n_states - skip_n_states);

        for (int i = 0; i < N; ++i) {
            double d = ie(i) - value_base;
            double d_abs = abs(d);
            if (d_abs > 0.005)
                num_above_error++;
            if (d_abs > largest_distance) {
                largest_distance = d_abs;
                largest_distance_i = i;
            }
            if (d_abs < abs(min_distance[i]))
                min_distance[i] = d;
            distance[i] = d;
        }

        if (largest_distance < best_distance) {
            best_distance = largest_distance;
            best_Ji = J_i;
        }

        std::cout << string_format(
                "G: %3.2f Step: %3d Num above: %3d Node: %2d Largest distance: %1.6f dist %7.6f delta: %.6f slow factor %f",
                model->get_param_value("G")[0],
                step,
                num_above_error,
                largest_distance_i,
                largest_distance,
                distance[largest_distance_i],
                delta[largest_distance_i],
                slow_factor)
                  << std::endl;

        if (num_above_error == 0) {
            break;
        }

        if (step > 0 && largest_distance > prev_largest_distance * 2.0)
            slow_factor *= 0.5;
        prev_largest_distance = largest_distance;

        for (int i = 0; i < N; ++i) {
            auto d_abs = abs(distance[i]);
            // auto delta_i = pow(d_abs, speed); // 0.003 * abs(d + 0.026) / 0.026
            float delta_i = slow_factor * d_abs;
            if (delta_i < 0.005)
                delta_i = 0.005;
            delta[i] = distance[i] > 0.0 ? delta_i : -delta_i;
        }

        J_i += delta;
        delete sim_result;
    }

    return {step, best_distance, best_Ji};
}


TArray1d eval_distance(const TArray1d &J_i, SimConfig &sim_config, int voi, float value_base) {
    sim_config.model()->set_param("J_i", J_i);
    auto [converged, sim_result] = tvb::simulate(sim_config, 1.0, voi);

    int N = sim_config.connectivity()->weights().rows();

    if (!converged) return TArray1d::Constant(N, 1e6);

    unsigned n_states = sim_result->getRecords().size();
    TArray1d distance(N);

    TArray1d total_ie = TArray1d::Zero(N);
    int skip_n_states = int(0.5 * n_states);
    for (unsigned j = skip_n_states; j < n_states; ++j)
        total_ie += sim_result->getRecords()[j].record.col(0);
    TArray1d ie = total_ie / (n_states - skip_n_states);

    delete sim_result;
    return ie - value_base;

}

class OptData {
public:
    OptData(const TArray1d &baseJi, Model *model, SimConfig *simConfig, int voi, float valueBase) : base_Ji(baseJi),
                                                                                                    model(model),
                                                                                                    sim_config(
                                                                                                            simConfig),
                                                                                                    voi(voi),
                                                                                                    value_base(
                                                                                                            valueBase) {}

    TArray1d base_Ji;
    Model *model;
    SimConfig *sim_config;
    int voi;
    float value_base;
};


double
sim_fn(const Eigen::VectorXd &vals_inp, Eigen::VectorXd *grad_out, void *opt_data) {
    OptData *odata = (OptData*)opt_data;

    TArray1d J_i = vals_inp(0) + vals_inp(1) * *(TArray1d *) opt_data;
    odata->model->set_param("J_i", J_i);
    auto [converged, sim_result] = tvb::simulate(*odata->sim_config, 1.0, odata->voi);

    int N = odata->sim_config->connectivity()->weights().rows();

    unsigned num_above_error = 0;
    unsigned n_states = sim_result->getRecords().size();
    double largest_distance = 0.0;
    unsigned largest_distance_i = 0;
    TArray1d distance(N);

    TArray1d total_ie = TArray1d::Zero(N);
    int skip_n_states = int(0.8 * n_states);
    for (unsigned j = skip_n_states; j < n_states; ++j)
        total_ie += sim_result->getRecords()[j].record.col(0);
    TArray1d ie = total_ie / (n_states - skip_n_states);

    for (int i = 0; i < N; ++i) {
        double d = ie(i) - odata->value_base;
        double d_abs = abs(d);
        if (d_abs > 0.005)
            num_above_error++;
        if (d_abs > largest_distance) {
            largest_distance = d_abs;
            largest_distance_i = i;
        }
        distance[i] = d;
    }

    float d = sqrt(distance.pow(2.0).sum()/N);
    std::cout << string_format("Evaluating a = %f, b = %f, distance = %f\n", vals_inp(0), vals_inp(1), d);

    return d;

}

//std::tuple<int, double, TArray1d> optimize_fic_Herzog_optimlib(SimConfig &sim_config, int voi, float value_base) {
//
//    int N = sim_config.connectivity()->weights().rows();
//
//    auto *model = sim_config.model();
//    TArray1d base_Ji = sim_config.connectivity()->weights().rowwise().sum();
//    float G = model->get_param_value("G")[0];
//    float avg_Ji = base_Ji.sum() / base_Ji.size();
//    base_Ji *= G / (4.0 * avg_Ji);
//
//    OptData opt_data{base_Ji, model, &sim_config, voi, value_base};
//
//    optim::ColVec_t initial(2);
//    initial << 1.0, 1.0;
//    bool success = optim::de(initial, sim_fn, &opt_data);
//
//    std::cout << success;
//
//}


std::tuple<bool, float, float, double, TArray1d> optimize_fic_Herzog(SimConfig &sim_config, int voi, float value_base, float a, float b) {
    int N = sim_config.connectivity()->weights().rows();

    auto *model = sim_config.model();
    TArray1d base_Ji = sim_config.connectivity()->weights().rowwise().sum();
    float G = model->get_param_value("G")[0];
    float avg_Ji = base_Ji.sum() / base_Ji.size();
    base_Ji -= avg_Ji;
    base_Ji *= G / (4.0 * avg_Ji);

    float da = 0.1;
    float db = 0.1;
    bool a_dir = true;
    float a_best, b_best;

    State initial_state;
    double prev_largest_distance = 0;
    // double sf = 1.0;
    double speed = 1.0;

    TArray1d J_i = model->get_param_value("J_i");
    TArray1d prev_J_i = J_i;
    float best_distance = 1e6;
    float current_best_distance = 1e6;
    TArray1d best_Ji;
    OptResult result;
    Monitor* sim_result;

    double slow_factor = 1.0;

    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<float> distr_a(1, 1.8);
    std::uniform_real_distribution<float> distr_b(0.3, 0.7);

    unsigned step = 0;
    int dir = 0;
    int n_fails = 0;
    int n_subs = 0;
    while (step < 1000) {
        float a0, b0;
        if (dir % 4 == 0 || dir % 4 == 2) {
            if (dir % 4 == 0)
                a0 = a + da;
            else
                a0 = a - da;
            b0 = b;
        } else {
            if (dir % 4 == 1)
                b0 = b + db;
            else
                b0 = b - db;
            a0 = a;
        }

//        a0 = distr_a(generator);
//        b0 = distr_b(generator);
        TArray1d d0 = eval_distance(a0 + b0 * base_Ji, sim_config, voi, value_base);

        float d0_f = d0.abs().maxCoeff();

        std::cout
                << string_format("G = %f (%f) at a=%f, b = %f, with %i evaluations from a=%f, b=%f (da=%.5f, db=%.5f)\n", G, d0_f, a0, b0,
                                 step, a, b, da, db);
        if (d0_f < current_best_distance) {
            current_best_distance = d0_f;
            if (current_best_distance < best_distance) {
                best_distance = current_best_distance;
                std::cout << string_format("G = %f NEW MINIMUM = %f at a=%f, b = %f\n", G, d0_f, a0, b0);
                a_best = a0;
                b_best = b0;
            }
            da = da < 0.1 ? da * 2.0 : 0.1;
            db = db < 0.1 ? db * 2.0 : 0.1;
            a = a0;
            b = b0;
            n_fails = 0;
            n_subs = 0;
        } else {
            n_fails++;
            if (n_subs < 6 && n_fails == 4) {
                da /= 2.0;
                db /= 2.0;
                n_subs++;
                n_fails = 0;
            } else {
                if (n_fails == 4 && n_subs == 6) {
                    n_fails = 0;
                    n_subs = 0;
                    a = distr_a(generator);
                    b = distr_b(generator);
                    da = 0.1;
                    db = 0.1;
                    current_best_distance = 1e6;
                    std::cout << string_format("G = %f NEW START at a=%f, b = %f\n", G, a, b);
                }
            }
            dir++;
        }

        delete sim_result;

        step += 1;

        if (best_distance < 0.005)
            break;
    }

    return {step < 1000, a_best, b_best, best_distance, a_best + b_best*base_Ji};
}

