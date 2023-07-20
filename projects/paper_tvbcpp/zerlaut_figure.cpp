//
// Created by imartin on 25-Oct-22.
//

#include <chrono>

#include <tvb-cpp/simulator/models/zerlaut.h>
#include <tvb-cpp/matplotlibcpp.h>
#include <tvb-cpp/simulator/stimulus.h>
#include <tvb-cpp/tools/csv_tools.h>
#include <tvb-cpp/datatypes/connectivity.h>
#include <tvb-cpp/simulator/integrators/euler_deterministic.h>
#include <tvb-cpp/simulator/coupling.h>
#include <tvb-cpp/simulator/simulate.h>
#include <tvb-cpp/simulator/integrators/euler_stochastic.h>
#include <tvb-cpp/simulator/modifier.h>
#include <tvb-cpp/tools/json.h>

namespace plt = matplotlibcpp;
using namespace tvb;

class ExpModifier : public tvb::ParamModifier {
    Float m_a, m_b, m_c, m_s;
    int m_n;
public:
    ExpModifier(int n, Float a, Float b, Float c, Float s): tvb::ParamModifier(n) {
        m_a = a;
        m_b = b;
        m_c = c;
        m_s = s;
    }

    TArray1d operator()(float t) override {
        return f(t);
    }

    TArray1d f(Float t) {
        Float r;
        if (t < m_a || t > m_c) r = 0.0;
        else if (t < m_b) {
            Float a = m_b - m_a;
            Float tt = 2.0 * (t-m_b) / a + 1.0;
            r = 1.0 + sin(tt * M_PI_2);
        }
        else {
            Float b = m_c - m_b;
            Float tt = 2.0 * (t-m_b) / b + 1.0;
            r = 1.0 + sin(tt * M_PI_2);
        }
        return TArray1d::Constant(m_n, 0.0025); //  + r/m_s);
    }
};


void sim_single() {
    Float simTime = 5000.0;
    Float dt = 0.1;

    // auto *model = new tvb::ZerlautAdaptationFirstOrder(1);
    auto *model = new tvb::ZerlautAdaptationSecondOrder(1);
    model->configure();
    model->set_param("K_ext_i", 100.0);
    model->set_param("T", 5.0);
    model->set_param("external_input_ex_ex", 2.5e-3);
    model->set_param("external_input_ex_in", 0.0);
    model->set_param("external_input_in_ex", 2.5e-3);
    model->set_param("external_input_in_in", 0.0);
    model->set_param("weight_noise", 1e-4);
    model->set_param("tau_OU", 5.0);

    TArray1d P_e(10);
    P_e << -0.0498,
            0.00506,
            -0.025,
            0.0014,
            -0.00041,
            0.0105,
            -0.036,
            0.0074,
            0.0012,
            -0.0407;
    TArray1d P_i(10);
    P_i << -0.0514,
            0.004,
            -0.0083,
            0.0002,
            -0.0005,
            0.0014,
            -0.0146,
            0.0045,
            0.0028,
            -0.0153;

    model->set_param("P_e", P_e);
    model->set_param("P_i", P_i);
    model->set_param("b_e", 60.0);


    model->init_dependant();
    for (auto p: model->get_param_list())
        std::cout << p << ": " << model->get_param_value(p)[0] << std::endl;

    auto mod = new ExpModifier(1, 375, 400, 750, 100.0);
    auto *modmod = new ModelModifier();
    modmod->configure(0.0, simTime, dt);
    modmod->addModifier("external_input_ex_ex", mod);
    // modmod->addModifier("external_input_ex_in", mod);

    std::vector<tvb::TArray2d> records;
    std::vector<tvb::TArray1d> vst;
    tvb::TArray2d state = tvb::TArray2d::Zero(1, model->n_vars());
    // state(0, 0) = 4.0;
    int n_steps = simTime / dt;
    for (int step = 0; step < n_steps; ++step) {
        auto d_state = (*model)(state, tvb::TArray2d::Zero(1, model->n_vars()), tvb::TArray1d::Zero(1));
        // d_state(7) += 1;
        state += dt * d_state;
        if (step%10 == 0) {
            records.push_back(state);
            vst.push_back((*mod)(step*dt));
        }
        // modmod->update(step, model);
    }

    int n_records = records.size();

    int t0 = 500;
    std::vector<std::vector<tvb::Float>> y_plot(model->n_vars());
    for (int n = 0; n < model->n_vars(); ++n)
        for (int step = t0; step < n_records; ++step) {
            y_plot[n].push_back(records[step](0, n));
        }

    std::vector<tvb::Float> s_plot;
    for (int step = t0; step < n_records; ++step) {
        s_plot.push_back(vst[step][0]);
    }

    auto ls = tvb::range(tvb::Float(t0*dt), simTime, n_records-1-t0);

    plt::figure_size(800, 1000);
    for (int n = 0; n < model->n_vars(); ++n) {
        plt::subplot(4, 2, n+1);
        plt::plot(ls, y_plot[n], {{"label", model->state_vars()[n]}});
        plt::legend();
        // plt::tight_layout();
    }

    // plt::plot(ls, s_plot, {{"color", "black"}, {"linestyle", "--"}});

    // Plot a red dashed line from given x and y data.
    // plt::plot(x, w,"r--");
    // Plot a line whose name will show up as "log(x)" in the legend.

//    plt::title("Zerlaut firing rate model");
//    plt::ylabel("KHz");
//    plt::xlabel("ms");
    // plt::ylim(0, 150);
    // Save the image (file format is determined by the extension)
    plt::save("./zerlaut_state_vars.png", 300);

}

double n_pdf(double x, double mean, double sigma) {
    return (1.0/(sigma*sqrt(2.0*M_PI))) * exp(-0.5*pow((x-mean)/sigma, 2.0));
}

void sim_fig1() {
    Float simTime = 5000.0;
    Float dt = 0.1;

    auto bes = tvb::range(0.0, 150.0, 10); // 76);
    std::vector<Float> Es;
    std::vector<Float> Is;
    std::vector<Float> C_ee;
    std::vector<Float> C_ii;

    auto *model = new tvb::ZerlautAdaptationSecondOrder(1);
    model->configure();

    model->set_param("K_ext_i", 100.0);
    model->set_param("T", 5.0);
    model->set_param("external_input_ex_ex", 2.5e-3);
    model->set_param("external_input_ex_in", 0.0);
    model->set_param("external_input_in_ex", 2.5e-3);
    model->set_param("external_input_in_in", 0.0);
    model->set_param("weight_noise", 1e-4);
    model->set_param("tau_OU", 5.0);

    TArray1d P_e(10);
    P_e << -0.0498,
            0.00506,
            -0.025,
            0.0014,
            -0.00041,
            0.0105,
            -0.036,
            0.0074,
            0.0012,
            -0.0407;
    TArray1d P_i(10);
    P_i << -0.0514,
            0.004,
            -0.0083,
            0.0002,
            -0.0005,
            0.0014,
            -0.0146,
            0.0045,
            0.0028,
            -0.0153;

    model->set_param("P_e", P_e);
    model->set_param("P_i", P_i);

    model->init_dependant();


    for (float  be: bes) {
        // auto *model = new tvb::ZerlautAdaptationFirstOrder(1);

        std::cout << "Simulating for be = " << be << std::endl;

        model->set_param("b_e", be);
        model->init_dependant();

        std::vector<tvb::TArray2d> records;
        std::vector<tvb::TArray1d> vst;
        tvb::TArray2d state = tvb::TArray2d::Zero(1, model->n_vars());
        // state(0, 0) = 4.0;
        int n_steps = simTime / dt;
        for (int step = 0; step < n_steps; ++step) {
            auto d_state = (*model)(state, tvb::TArray2d::Zero(1, model->n_vars()), tvb::TArray1d::Zero(1));
            state += dt * d_state;
            if (step % 10 == 0) {
                records.push_back(state);
            }
        }
        Es.push_back(state(0, 0));
        Is.push_back(state(0, 1));
        C_ee.push_back(state(0, 2));
        C_ii.push_back(state(0, 4));
    }

    plt::plot(bes, Es, {{"color", "green"}, {"label", "excitatory"}});
    plt::plot(bes, Is, {{"color", "red"}, {"label", "inhibitory"}});
    plt::legend();

    // Plot a red dashed line from given x and y data.
    // plt::plot(x, w,"r--");
    // Plot a line whose name will show up as "log(x)" in the legend.

//    plt::title("Zerlaut firing rate model");
//    plt::ylabel("KHz");
//    plt::xlabel("ms");
    // plt::ylim(0, 150);
    // Save the image (file format is determined by the extension)
    plt::save("./zerlaut_fig1.png", 300);


    plt::clf();
    auto mu = Es[4];
    auto sigma = sqrt(C_ee[4]);
    auto x = tvb::range(mu - 3*sigma, mu + 3*sigma, 100);
    std::vector<double> y(x.size());
    std::transform(x.begin(), x.end(), y.begin(), [mu, sigma](double x) { return n_pdf(x, mu, sigma); });
    plt::plot(x, y, {{"color", "green"}, {"label", "excitatory"}});
    mu = Is[4];
    sigma = sqrt(C_ii[4]);
    x = tvb::range(mu - 3*sigma, mu + 3*sigma, 100);
    std::transform(x.begin(), x.end(), y.begin(), [mu, sigma](double x) { return n_pdf(x, mu, sigma); });
    plt::plot(x, y, {{"color", "red"}, {"label", "inhibitory"}});
    plt::legend();
    plt::save("./zerlaut_fig1_b.png", 300);

}

void sim_whole() {
    tvb::TArray2d C = tvb::csv_load("0001_1_Counts.csv");

    int N = C.rows();

    C = C / C.rowwise().sum().maxCoeff() * 2.0;
    // tvb::csv_save("sc_d_norm.csv", C);

    tvb::TArray2d tl = tvb::TArray2d::Zero(C.rows(), C.cols());
    tvb::Connectivity con(C, tl, 1e6);

    std::chrono::milliseconds total_time(0);
    std::cout << string_format("Starting computation for: %s", "C") << std::endl;

    //auto *model = new tvb::Montbrio(N, rp.t_start, rp.t_end, rp.dt);
    auto *model = new tvb::ZerlautAdaptationSecondOrder(N);

    float dt = 0.1;

    tvb::TArray1d sigmas(8);
    sigmas << 1e-6, 1e-6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    auto *integrator = new tvb::EulerStochastic(dt, new tvb::Additive(sigmas, 0.1));
    // auto *integrator = new tvb::EulerDeterministic();
    auto *coupling = new tvb::CouplingLinearSparse(con.weights(), con.delays(), model->cvars());

    auto start = std::chrono::high_resolution_clock::now();
    SimConfig sim_config;

    sim_config.setModel(model);
    sim_config.setConnectivity(&con);
    sim_config.setIntegrator(integrator);
    sim_config.setCoupling(coupling);
    sim_config.setIntegrationInterval(0, 3000);
    sim_config.setNumIterations(1);
    sim_config.setTimeDelta(dt);
    sim_config.setDeltaIntegration(0.00001);

    auto [converged, monitor] = tvb::simulate(sim_config, 1.0, 0);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start);

    std::cout << string_format("Simulation time (%s): %d msecs", "C", duration.count()) << std::endl;

    size_t n_records = monitor->getRecords().size();
    std::vector<std::vector<Float>> y_plot(N, std::vector<Float>(n_records));
    for (unsigned t = 0; t < n_records; ++t)
        for (unsigned n = 0; n < N; ++n)
            y_plot[n][t] = 1e-3 * monitor->getRecords()[t].record(n, 0);

    // tvb::csv_save("./paper_RWW_BOLD_TVBCPP.csv", y_plot);

    // Plot line from given x and y data. Color is selected automatically.
    std::vector<Float> ls(n_records);
    std::transform(monitor->getRecords().begin(), monitor->getRecords().end(), ls.begin(),
                   [](const Monitor::Record &r) { return r.time/1000; });

    for (int i = 0; i < N; ++i)
        plt::plot(ls, y_plot[i]);

    // Plot a red dashed line from given x and y data.
    // plt::plot(x, w,"r--");
    // Plot a line whose name will show up as "log(x)" in the legend.

    plt::title("Zerlaut whole brain");
//    plt::ylim(0.0, 30.0);
    plt::ylabel("ve (Hz)");
    plt::xlabel("Seconds");
    // Save the image (file format is determined by the extension)
    plt::save("./zerlaut_whole.png", 300);
}


int main(int argc, char ** argv) {
    // sim_single();
    // sim_whole();
    sim_fig1();
}
