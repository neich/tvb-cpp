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

#ifndef TVB_CPP_SIMULATE_FCD_H
#define TVB_CPP_SIMULATE_FCD_H

#include <definitions.h>
#include <simulator/simulate.h>

class SimulateFCD {
protected:
    double dtt = 1e-3; // Sampling rate of simulated neuronal activity (seconds)
    double dt = 0.1;

    double TR = 2.0;
    double Tmax = 220.0;
    double Toffset = 10.0;
    double Tmaxneuronal;

    int windowSize = 10;
    int N_windows;

    void recomputeTmaxneuronal() {
        Tmaxneuronal = int((Tmax+Toffset)*(TR/dtt));
        N_windows = int(round((Tmax-windowSize)/3.0));
    }

public:
    SimulateFCD() {
        recomputeTmaxneuronal();
    }

    SimulateFCD(double dtt, double dt,
                double tr, double tmax, double toffset,
                int windowSize) : dtt(dtt), dt(dt),
                                  TR(tr), Tmax(tmax),
                                  Toffset(toffset),
                                  windowSize(windowSize) {
        recomputeTmaxneuronal();
    }

    tvb::Matrixd computeSubjectBold(const tvb::StateTrack& signal, double dt, const tvb::Vectori& areasToSimulate={}) const;
    tvb::Matrixd simulateSingleSubject(tvb::SimConfig &simConfig) const;

};


#endif //TVB_CPP_SIMULATE_FCD_H
