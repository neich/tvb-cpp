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

#include <chrono>

#include <fic/functions/g_optim.h>
#include <simulator/models/reduced_ww_ext.h>
#include <fic/functions/balance_fic.h>

using namespace std;
using namespace tvb;

TArray2dMap processBOLDSignals(const vector<tvb::TArray2d>& BOLDsignals,
                               const DistanceSettings& distanceSettings,
                               const Filter& filter) {

    int NumSubjects = BOLDsignals.size();
    int N = BOLDsignals[0].rows(); // get the first key to retrieve the value of N = number of areas

    // # First, let's create a data structure for the distance measurement operations...
    TArray2dMap measureValues;
    for (auto const& [ key, measure ] : distanceSettings) {
        measure.init(NumSubjects, N);
    }

    // # Loop over subjects
    int nsub = 0;
    for (const auto & signal : BOLDsignals) {
        // if verbose: print('   BOLD {}/{} Subject: {} ({}x{})'.format(pos, NumSubjects, s, BOLDsignals[s].shape[0], BOLDsignals[s].shape[1]), end='', flush=True)
        // signal = BOLDsignals[s]  # LR_version_symm(tc[s])

        cout << "Processing subject " << nsub << endl;
        for (auto &[ds, measure] : distanceSettings) { // # Now, let's compute each measure and store the results
            TArray2d procSignal = measure.from_fMRI(signal);
            measure.accumulate(procSignal, nsub);
        }
        nsub++;
    }

    for (auto const& [ds, measure] : distanceSettings) { //  # finish computing each distance measure
        measureValues[ds] = measure.postprocess();
    }

    return measureValues;
}

TArray2dMap distanceForOne_G(double we, const tvb::TArray1d& J_i,
                             SimConfig &sim_config, int N, int NumSimSubjects,
                             const SimulateFCD &sim_fcd,
                             const DistanceSettings &distanceSettings) {

    ReducedWongWangExcInh *model = dynamic_cast<ReducedWongWangExcInh*>(sim_config.model());
    model->G.fill(we);
    model->J_i = J_i;

    cout << string_format("   --- BEGIN TIME @ we=%f ---", we) << endl;
    auto start = std::chrono::high_resolution_clock::now();

    vector<tvb::TArray2d> simulatedBOLDs;
    for (int nsub = 0; nsub < NumSimSubjects; ++nsub) { //  # trials. Originally it was 20.
        cout << string_format("   Simulating we=%f -> subject %d/%d!!!", we, nsub, NumSimSubjects) << endl;
        TArray2d bold_signal = sim_fcd.simulateSingleSubject(sim_config);
//        tvb::Bold bold_monitor(simConfig, {3});
//        StateTrack bold_result = bold_monitor.apply(result.m_times, result.m_states);
//        TArray2d bold_signal = stateTrackToMatrix(bold_result);
        simulatedBOLDs.push_back(bold_signal);
    }

    TArray2dMap dist = processBOLDSignals(simulatedBOLDs, distanceSettings);
    // dist["We"] = we;
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    cout << string_format("   --- TOTAL TIME: %d seconds ---\\n", duration);
    return dist;
}
