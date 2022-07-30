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

#include <simulator/simulator.h>
#include "simulator/integrators/euler_deterministic.h"
#include "simulator/models/zerlaut.h"
#include "simulator/models/reduced_ww_ext.h"
#include <matplotlibcpp.h>
#include <chrono>

using namespace tvb;
namespace plt = matplotlibcpp;

int main(int /* argc */, char ** /* argv */ ) {

    auto *integrator = new tvb::EulerDeterministic();

    double start_time = 0;
    double end_time = 10000;
    double dt = 0.1;
    integrator->configure(start_time, end_time, dt);

    auto n_steps = int((end_time - start_time) / dt);
    auto n_reg = 1;

    // Model *model = new ZerlautFirstOrder(n_reg);
    // ZerlautAdptationSecondOrder *model = new ZerlautAdptationSecondOrder(n_reg);
//    model->external_input.fill(0.1);
    auto *model = new ReducedWongWangExcInh(n_reg);

    State state;
    model->initial(state);
    StateTrack *result = model->create_track();
    double t = start_time;

    // result->push(state, t);

    TArray1d node_coupling = TArray1d::Zero(n_reg);
    TArray1d local_coupling = TArray1d::Zero(n_reg, model->n_vars());
    TArray1d stimulus = TArray1d::Constant(model->n_vars(), 0);

    auto start = std::chrono::high_resolution_clock::now();
    for (int step = 1; step <= n_steps; ++step) {
        state = integrator->scheme(state, *model, node_coupling, local_coupling, stimulus);
        t += integrator->dt();
        result->push(state, t);
    }
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
    std::cout << string_format("Simulation time: %d msecs", duration.count());

    int N = result->m_states[0].rows();

    size_t t_max = result->m_states.size();

    std::vector<std::vector<double>> y_plot(model->n_vars(), std::vector<double>(t_max));
    for (unsigned t = 0; t < t_max; ++t)
        for (unsigned r = 0; r < model->n_vars(); ++r)
            y_plot[r][t] = result->m_states[t](0, r);

    plt::figure_size(1200, 780);
    // Plot line from given x and y data. Color is selected automatically.
    plt::plot(result->m_times, y_plot[0]);
    plt::plot(result->m_times, y_plot[1]);
    // Plot a red dashed line from given x and y data.
    // plt::plot(x, w,"r--");
    // Plot a line whose name will show up as "log(x)" in the legend.

    plt::title("Sample figure");
    // Enable legend.
    plt::legend();
    // Save the image (file format is determined by the extension)
    plt::save("./test_model.png");
}


