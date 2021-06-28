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

#include <simulator/simulate.h>
#include <simulator/models/reduced_ww_ext.h>

#include <fic/functions/balance_fic.h>
#include <matplotlibcpp.h>

using namespace tvb;
namespace plt = matplotlibcpp;

OptResult optimize_fic(double G, SimConfig &sim_config) {

    int N = sim_config.connectivity()->weights().rows();

    ReducedWongWangExcInh *model = dynamic_cast<ReducedWongWangExcInh *>(sim_config.model());
    model->G.fill(G);

    double be_ae = model->b_e(0) / model->a_e(0);

    Vectord delta(N);
    delta.fill(1.0);
    Vectord distance(N);
    distance.fill(10.0);
    State initial_state;
    double prev_largest_distance = 0;
    // double sf = 1.0;
    double speed = 1.0;

    Vectord J_i = model->J_i;
    Vectord prev_J_i = model->J_i;
    // double best_distance = 1e6;
    Vectord best_Ji;
    OptResult result;
    StateTrack sim_result;

    double slow_factor = 1.0;

    for (unsigned step = 0; step < 100; ++step) {
    //for (unsigned step = 0; step < 500; ++step) {
        model->J_i = J_i;
        sim_result = tvb::simulate(sim_config);

//        size_t t_max = sim_result.m_states.size();
//
//        std::vector<std::vector<double>> y_plot(N, std::vector<double>(t_max));
//        for (unsigned t = 0; t < t_max; ++t)
//            for (unsigned r = 0; r < N; ++r)
//                y_plot[r][t] = sim_result.m_states[t](r, 3);
//
//        plt::figure_size(1200, 780);
//        // Plot line from given x and y data. Color is selected automatically.
//        for (unsigned i = 0; i < N; ++i)
//            plt::plot(sim_result.m_times, y_plot[i]);
//        // Plot a red dashed line from given x and y data.
//        // plt::plot(x, w,"r--");
//        // Plot a line whose name will show up as "log(x)" in the legend.
//
//        plt::title("Sample figure");
//        // Enable legend.
//        plt::legend();
//        plt::ylim(0.0, 0.5);
//        // Save the image (file format is determined by the extension)
//        plt::save("./sim_test.png");


        unsigned num_above_error = 0;
        unsigned n_states = sim_result.m_states.size();
        double largest_distance = 0.0;
        unsigned largest_distance_i = 0;
        Vectord distance(N);
        Vectord min_distance(N);

        Vectord total_ie = Vectord::Zero(N);
        int skip_n_states = int(0.2 * n_states);
        for (unsigned j = skip_n_states; j < n_states; ++j)
            total_ie += sim_result.m_states[j].col(3);
        Vectord ie = total_ie / (n_states - skip_n_states);

        for (int i = 0; i < N; ++i) {
            double d = ie(i) - be_ae + 0.026;
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

        std::cout << string_format(
                "G: %3.2f Step: %3d Num above: %3d Node: %2d Largest distance: %1.6f dist %7.6f delta: %.6f speed %f",
                G,
                step,
                num_above_error,
                largest_distance_i,
                largest_distance,
                distance[largest_distance_i],
                delta[largest_distance_i],
                speed)
                  << std::endl;

        if (num_above_error == 0) {
            break;
        }

        if (step > 0 && largest_distance > prev_largest_distance * 2.0)
            speed *= 0.9;
        prev_largest_distance = largest_distance;

        for (int i = 0; i < N; ++i) {
            auto d_abs = abs(distance[i]);
            auto delta_i = pow(d_abs, speed); // 0.003 * abs(d + 0.026) / 0.026
            delta_i = slow_factor * d_abs;
            delta[i] = distance[i] > 0.0 ? delta_i : -delta_i;
        }

//        for (int i = 0; i < N; ++i) {
//            auto d_abs = abs(distance[i]);
//            if (d_abs > 0.005) {
//                auto delta_i = slow_factor * d_abs / 0.1;
//                if (delta_i < 0.005)
//                    delta_i = 0.0;
//            }
//            delta[i] = distance[i] > 0.0 ? delta_i : -delta_i;
//        }


        J_i += delta;

    }

    result.m_states = sim_result.m_states;
    result.m_times = sim_result.m_times;
    result.m_Jis = J_i;

    return result;
}

