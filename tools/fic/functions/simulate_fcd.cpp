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

#include <simulator/monitors/bold_Stephan2007.h>

#include <external/numpy/numpy.h>
#include <fic/functions/simulate_fcd.h>

using namespace tvb;

Matrixd SimulateFCD::computeSubjectBold(const StateTrack& signal, double dt, const Vectori& areasToSimulate) const {
    BoldStephan2007 bold(signal.m_states.size()*dtt,signal.m_states[0].rows(), dtt, {2});
    Matrixd result = bold.apply(stateTrackToMatrix(signal, 2));
    int step = int(round(TR/dtt));
    return result(Eigen::all, Eigen::seq(step-1, Eigen::last, step));
}


Matrixd SimulateFCD::simulateSingleSubject(SimConfig &simConfig) const{
    simConfig.setIntegrationInterval(0.0, Tmaxneuronal);
    simConfig.setTimeDelta(dt);
    simConfig.setSamplingRate(10);
    StateTrack sresult = simulate(simConfig);
    Matrixd bds = computeSubjectBold(sresult, simConfig.dt());
    return bds;
}

