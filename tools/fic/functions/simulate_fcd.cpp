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
#include "simulator/monitors/bold_tvb.h"

using namespace tvb;


TArray2d SimulateFCD::simulateSingleSubject(SimConfig &simConfig, int voi) const{
    simConfig.setIntegrationInterval(0.0, Tmaxneuronal);
    simConfig.setTimeDelta(dt);
    simConfig.setSamplingRate(10);
    int  n_roi = simConfig.connectivity()->weights().rows();
    auto *monitor = new BoldTVB(n_roi, TR, dt, {voi});
    simulate(simConfig, 1.0, voi);
    TArray2d bds(n_roi, monitor->getRecords().size());
    int n = 0;
    for (auto const& r: monitor->getRecords())
        bds.col(n++) = r.record.col(0);
    return bds;
}

