//
// Created by imartin on 25-Oct-22.
//

#include <chrono>

#include <tvb-cpp/simulator/models/zerlaut.h>
#include <tvb-cpp/matplotlibcpp.h>
#include <tvb-cpp/simulator/stimulus.h>
#include <tvb-cpp/tools/csv_tools.h>
#include <tvb-cpp/tools/npz_tools.h>
#include <tvb-cpp/datatypes/connectivity.h>
#include <tvb-cpp/simulator/integrators/euler_deterministic.h>
#include <tvb-cpp/simulator/coupling.h>
#include <tvb-cpp/simulator/simulate.h>
#include <tvb-cpp/simulator/integrators/euler_stochastic.h>
#include <tvb-cpp/simulator/modifier.h>
#include <tvb-cpp/tools/json.h>
#include <tvb-cpp/tools/threadpool.h>

namespace plt = matplotlibcpp;
using namespace tvb;

void configure_model_divolo(ZerlautAdaptationSecondOrder* model) {
    model->configure();
    model->set_param("C_m", 200);
    model->set_param("K_ext_i", 100.0);
    model->set_param("T", 5.0);
    model->set_param("weight_noise", 1e-4);
    model->set_param("tau_OU", 5.0);
    model->set_param("b_e", 60.0);

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
    model->set_param("external_input_ex_ex", 2.5e-3);
    model->set_param("external_input_ex_in", 0.0);
    model->set_param("external_input_in_ex", 2.5e-3);
    model->set_param("external_input_in_in", 0.0);


    model->init_dependant();

}

RunParams run(RunParams rp) {

    string f_prefix = (path(rp.path_out) / path(rp.file_prefix)).lexically_normal().string();
    for (auto const &p: rp.params) {
        f_prefix += string_format("_%s_%.2f", p.name.c_str(), p.value);
    }
    string filename = f_prefix + ".png";

    if (!rp.force_output && std::filesystem::exists(filename)) {
        std::cout << string_format("File %s already exists", filename.c_str()) << std::endl;
        delete rp.monitor;
        rp.monitor = nullptr;
        return rp;
    }

    tvb::TArray2d C;
    if (rp.file_weights.ends_with(".csv"))
        C = tvb::csv_load(rp.file_weights);
    else if (rp.file_weights.ends_with(".npz"))
        C = tvb::npz2Matrixd(rp.file_weights, "SC");
    else
        throw std::runtime_error(string_format("Unknown file extension for: %s", rp.file_weights.c_str()));

    int N = C.rows();

    Float k = 0.15 / (C.rowwise().sum().sum() / N);
    // C *= k;
    // tvb::csv_save("sc_d_norm.csv", C);

    tvb::TArray2d tl;
    if (!rp.file_lengths.empty())
        tl = tvb::csv_load(rp.file_lengths);
    else
        tl = tvb::TArray2d::Zero(C.rows(), C.cols());
    auto *con = new tvb::Connectivity(C, tl, rp.speed);

    milliseconds total_time(0);
    std::cout << string_format("Starting computation for: %s", filename.c_str()) << std::endl;

    // auto *model = new tvb::Montbrio(N, rp.t_start, rp.t_end, rp.dt);
    // auto *model = new tvb::ReducedWongWangExcInh(N);
    auto *model = new ZerlautGABA(N);
    for (auto const &p: rp.params)
        if (std::isalpha(p.name[0])) model->set_param(p.name, p.value);

    // rp.monitor = new tvb::BoldTVB(N, 720.0, rp.dt, {0});
    // rp.monitor = new tvb::BoldBalloonWindkessel(N, 1.0, 720.0, rp.dt, {0});
    // rp.monitor = new tvb::RawSubSample(1.0, rp.dt, {3});

    TArray1d sigmas = TArray1d::Constant(model->n_vars(), 0.0);
    for (auto const &p: rp.params)
        if (p.name[0] == '_') {
            auto idx = std::stoi(p.name.substr(2, 1));
            sigmas[idx] = p.value;
        }


    // sigmas << 3e-5, 3e-5, 0.0, 0.0;
    auto *integrator = new tvb::EulerStochastic(rp.dt, new Additive(sigmas, rp.dt));
    // auto *integrator = new tvb::EulerDeterministic();

    auto coupling = new tvb::CouplingLinearSparse(con->weights(), con->delays(), model->cvars());

    auto start = std::chrono::high_resolution_clock::now();

    SimConfig sim_config;

    sim_config.setModel(model);
    sim_config.setConnectivity(con);
    sim_config.setIntegrator(integrator);
    // sim_config.setMonitor(rp.monitor);
    sim_config.setCoupling(coupling);
    sim_config.setIntegrationInterval(rp.t_start, rp.t_end);
    sim_config.setNumIterations(1);
    sim_config.setDeltaIntegration(0.00001);

    //auto [step, distance, J_i] = optimize_fic(sim_config, rp.voi, rp.value_base);
    auto [found, a, b, distance, J_i] = optimize_fic_Herzog(sim_config, rp.voi, rp.value_base);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start);

    std::cout << string_format("Optimization time (%s): %d msecs", filename.c_str(), duration.count()) << std::endl;

    total_time += duration;

//    std::vector<long unsigned> shape{J_i.rows(), J_i.cols()};
//    npy::SaveArrayAsNumpy(f_prefix + ".npz", false, shape.size(), shape.data(), (Float*)J_i.data());

    if (found) {

        TArray1d2npz(J_i, f_prefix + ".npz", "J_i");
        TArray1d2npz(TArray1d::Constant(1, 1, model->get_param_value("G")[0]), f_prefix + ".npz", "G");
        TArray1d2npz(sigmas, f_prefix + ".npz", "s");

        model->set_param("J_i", J_i);
        auto [converged, sim_result] = tvb::simulate(sim_config, 1.0, rp.voi);

        save_fig(sim_result, f_prefix);

        std::ofstream out_txt(f_prefix + ".txt");
        out_txt << string_format("MINIMUM = %f at a=%f, b = %f\n", distance, a, b);


        delete sim_result;
        rp.file_out = filename;
    }

    delete model;
    delete coupling;

    return rp;
}

void sim_whole() {
    tvb::TArray2d C = tvb::npy2Matrixd("/mnt/c/Users/IMARTIN/Dropbox/work/git/research/neuro/data/anesthesia/SampleData_Propofol/SC_Controls_CleanAverage_30perc.npy");

    int N = C.rows();

    C = C / (C.rowwise().sum().maxCoeff() * 100.0);
    // tvb::csv_save("sc_d_norm.csv", C);

    // C.setZero();

    tvb::TArray2d tl = tvb::TArray2d::Zero(C.rows(), C.cols());
    tvb::Connectivity con(C, tl, 1e6);

    std::chrono::milliseconds total_time(0);
    std::cout << string_format("Starting computation for: %s", "C") << std::endl;

    //auto *model = new tvb::Montbrio(N, rp.t_start, rp.t_end, rp.dt);
    auto *model = new tvb::ZerlautAdaptationSecondOrder(N);

    Float simTime = 5000.0;
    Float dt = 0.1;

    configure_model_divolo(model);


    tvb::TArray1d sigmas(8);
    sigmas << 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 1e-5;
    //auto *integrator = new tvb::EulerStochastic(dt, new tvb::Additive(sigmas, 0.1));
    auto *integrator = new tvb::EulerDeterministic(dt);
    auto *coupling = new tvb::CouplingLinearSparse(con.weights(), con.delays(), model->cvars());

    tvb::TArray1d Gs(8);
    sigmas << 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 9.0, 15;
    tvb::ThreadPool<RunParams> tp(sigmas.size());
    tp.start();

    for (auto G: Gs) {
        tp.queue_job([G] {
        // auto start = std::chrono::high_resolution_clock::now();
        SimConfig sim_config;

        sim_config.setModel(model);
        sim_config.setConnectivity(&con);
        sim_config.setIntegrator(integrator);
        sim_config.setCoupling(coupling);
        sim_config.setIntegrationInterval(0, simTime);
        sim_config.setNumIterations(1);
        sim_config.setTimeDelta(dt);
        sim_config.setDeltaIntegration(0.00001);

        auto [converged, monitor] = tvb::simulate(sim_config, 1.0, 0);

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - start);

        // std::cout << string_format("Simulation time (%s): %d msecs", "C", duration.count()) << std::endl;

        size_t n_records = monitor->getRecords().size();
        int skip = 0;
        int nplot = n_records - skip;
        std::vector<std::vector<Float>> y_plot(N, std::vector<Float>(nplot));
        for (unsigned t = skip; t < n_records; ++t)
            for (unsigned n = 0; n < N; ++n)
                y_plot[n][t - skip] = monitor->getRecords()[t].record(n, 0);

        // tvb::csv_save("./paper_RWW_BOLD_TVBCPP.csv", y_plot);

        // Plot line from given x and y data. Color is selected automatically.
        std::vector<Float> ls(nplot);
        for (int i = skip; i < n_records; ++i) {
            ls[i - skip] = monitor->getRecords()[i].time / 1000;
        }
//    std::transform(monitor->getRecords().begin(), monitor->getRecords().end(), ls.begin(),
//                   [](const Monitor::Record &r) { return r.time/1000; });

        for (int i = 0; i < N; ++i)
            plt::plot(ls, y_plot[i]);

        // Plot a red dashed line from given x and y data.
        // plt::plot(x, w,"r--");
        // Plot a line whose name will show up as "log(x)" in the legend.

        plt::title("Zerlaut whole brain");
        plt::ylim(0.0, 0.025);
        plt::ylabel("ve (Hz)");
        plt::xlabel("Seconds");
        // Save the image (file format is determined by the extension)
        plt::save("./zerlaut_whole.png", 300);
        });
    }
}


int main(int argc, char ** argv) {
    sim_whole();
}
