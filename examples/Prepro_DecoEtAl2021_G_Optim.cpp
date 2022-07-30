//// ==========================================================================
//// ==========================================================================
////  Computes the Functional Connectivity Dynamics (FCD)
////
////  From the original code:
//// --------------------------------------------------------------------------
////  OPTIMIZATION GAIN
////
////  Taken from the code (FCD_LSD_empirical.m) from:
////  [DecoEtAl_DRAFT] Gustavo Deco, Kevin Aquino, Aurina Arnatkeviciute, Stuart Oldham, Kristina Sabaroedin,
////  Nigel Rogasch, Morten L. Kringelbach, and Alex Fornito, "Transcriptomic heterogeneity gives rise to
////  cortical ignition and a hierarchy of timescales", DRAFT, to be submitted
////
////  Translated to Python & refactoring by Gustavo Patow
//// ==========================================================================
//// ==========================================================================

#include <string>
#include <unordered_map>

#include <definitions.h>
#include <tools/npz_tools.h>
#include <external/numpy/numpy.h>
#include <fic/functions/observers/sw_fcd.h>
#include <fic/functions/observers/gbc.h>
#include <fic/functions/g_optim.h>
#include <fic/functions/bold_filters.h>
#include <load_or_compute.h>
#include <simulator/models/reduced_ww_ext.h>
#include <fic/functions/balance_fic.h>
#include <simulator/integrators/euler_stochastic.h>


std::string baseInPath = "Data_Raw/DecoEtAl2020/";
std::string baseOutPath = "Data_Produced/DecoEtAl2020/";

using namespace std;
using namespace tvb;

void prepro() {

    int N = 68;
    int NSUB = 389;
    int NumTrials = 1; // 10

    // load genetic info
    cout << "Loading DKcortex_selectedGenes.npz" << endl;
    TArray2d expMeasures = npz2Matrixd(baseInPath + "DKcortex_selectedGenes.npz", "expMeasures");

    cout << "Rows: " << expMeasures.rows() << ", cols: " << expMeasures.cols() << endl;

    auto coefei = expMeasures(Eigen::all, Eigen::seq(17, 24)).rowwise().sum()
                  / expMeasures(Eigen::all, Eigen::seq(1, 5)).rowwise().sum();  // ampa+nmda/gaba
    auto ratioEI = TArray1d(N).setZero();
    ratioEI(Eigen::seq(0, coefei.size() - 1)) = coefei / (coefei.maxCoeff() - coefei.minCoeff());
    ratioEI = ratioEI - ratioEI.maxCoeff() + 1.0;
    ratioEI(Eigen::seq(35, 67)) = ratioEI(Eigen::seq(1, 33));

    cout << "loading SC_GenCog_PROB_30.npz" << endl;
    auto GrCV = npz2Matrixd(baseInPath + "SC_GenCog_PROB_30.npz", "GrCV");
    cout << "loading DKatlas_noGSR_timeseries.npz";
    std::vector<tvb::TArray2d> ts = npz2VecMatrixd(baseInPath + "DKatlas_noGSR_timeseries.npz", "ts");

    TArray2d C(N, N);
    C.topLeftCorner(34, 34) = GrCV(Eigen::seq(0, 33), Eigen::seq(0, 33));
    C.topRightCorner(34, 34) = GrCV(Eigen::seq(0, 33), Eigen::seq(41, 74));
    C.bottomLeftCorner(34, 34) = GrCV(Eigen::seq(41, 74), Eigen::seq(0, 33));
    C.bottomRightCorner(34, 34) = GrCV(Eigen::seq(41, 74), Eigen::seq(41, 74));
    C = 0.2 * C / C.maxCoeff();

    BandPassFilter bpf(0.008, 0.08, 0.754);
    SW_FC swFCD(80, 18, true, bpf);
    FunctionalConnectivityStandard FC;
    GBC_FC gbc;

    DistanceSettings distanceSettings = {{"FC",    FC},
                                         {"swFCD", swFCD},
                                         {"GBC",   gbc}};

    // Transform empirical subjects
    cout << "Transforming empirical subjects\n";
    TArray1di tcrange(N);
    tcrange << arange<int>(0, 34), arange<int>(41, 75);


    std::vector<TArray2d> tc_transf(NSUB);
    std::fill(tc_transf.begin(), tc_transf.end(), tvb::TArray2d::Zero(N, ts.size()));

    for (unsigned i = 0; i < NSUB; ++i)
        for (unsigned j = 0; j < ts.size(); ++j)
            tc_transf[i].col(j) = ts[j](tcrange, i);

    cout << "Processing empirical subjects\n";
    TArray2dMap FCemp_cotsampling =
            tvb::load_or_compute(baseOutPath + "fNeuro_emp.npz",
                                 [&tc_transf, &distanceSettings]() -> TArray2dMap {
                                     return processEmpiricalSubjects(tc_transf,
                                                                     distanceSettings);
                                 }
            );

    TArray2d FCemp = FCemp_cotsampling["FC"];
    TArray2d cotsampling = FCemp_cotsampling["swFCD"];
    TArray2d GBCemp = FCemp_cotsampling["GBC"];

    std::string J_fileNames = baseOutPath + "J_Balance_we%.1f.npz";
    std::string baseName = baseOutPath + "fitting_we%.1f.npz";

    // TArray1d WEs = arange<double>(0, 3, 0.1);
    TArray1d WEs = arange<double>(0, 0.001, 0.001); //  DEBUG only!!!
    int numWEs = WEs.size();
    TArray1d FCfitt = TArray1d::Zero(numWEs);
    TArray1d swFCDfitt = TArray1d::Zero(numWEs);
    TArray1d GBCfitt = TArray1d::Zero(numWEs);

    // Configure simulation
    tvb::SimConfig sim_config;

    tvb::TArray2d tl(N, N);
    tvb::generate(tl, []() { return (5.0 * (double) rand() / (RAND_MAX)); });
    tvb::Connectivity con(C, tl, 1e100);
    auto *model = new tvb::ReducedWongWangExcInh(N);
    // model->G.fill(we);
    TArray1d sigmas(4);
    sigmas << 3e-5, 3e-5, 0.0, 0.0;
    auto *integrator = new tvb::EulerStochastic(new Additive(sigmas, 0.1));

    sim_config.setModel(model);
    sim_config.setIntegrator(integrator);
    sim_config.setConnectivity(&con);
    sim_config.setHistory(new HistoryNoDelays(con.weights(), con.delays(), {3}));
    sim_config.setIntegrationInterval(0.0, 10000.0);
    sim_config.setTimeDelta(0.1);

    SimulateFCD sim_fcd(0.001, 0.1, 0.754, 616.0, 14.0, 10.0);

    for (unsigned pos = 0; pos < numWEs; ++pos) {
        double we = WEs[pos];


        TArray2d J_i =
                tvb::load_or_compute_index(string_format(J_fileNames, we), "J_i",
                                           [&we, &sim_config]() -> TArray2d {
                                               return optimize_fic(we, sim_config).m_Jis;
                                           });

        TArray2dMap FCsimul_cotsamplingsim =
                tvb::load_or_compute(string_format(baseName, we),
                                     [we, &sim_config, &J_i, &N, &NumTrials, &sim_fcd, &distanceSettings]() -> TArray2dMap {
                                         return distanceForOne_G(we, J_i,
                                                                 sim_config, N, NumTrials,
                                                                 sim_fcd,
                                                                 distanceSettings);
                                     });


        TArray2d FC_sim = FCsimul_cotsamplingsim["FC"];
        TArray2d swFCD_sim = FCsimul_cotsamplingsim["swFCD"];
        TArray2d GBC_sim = FCsimul_cotsamplingsim["GBC"];

        swFCDfitt[pos] = swFCD.distance(cotsampling, swFCD_sim);
        FCfitt[pos] = FC.distance(FCemp, FC_sim);
        GBCfitt[pos] = gbc.distance(GBCemp, GBC_sim);
        cout << "swFCDfitt = " << swFCDfitt[pos]
             << ", FCfitt = " << FCfitt[pos]
             << ", GBCfitt = " << FCfitt[pos] << endl;
    }
}

int main(int /* argc */, char ** /* argv */ ) {
    prepro();
}
