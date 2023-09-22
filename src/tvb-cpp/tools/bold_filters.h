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

#ifndef TVB_CPP_BOLD_FILTERS_H
#define TVB_CPP_BOLD_FILTERS_H

#include "tvb-cpp/definitions.h"

#include <tvb-cpp/tools/scipy/signal/filter_design.h>
#include <tvb-cpp/tools/scipy/signal/signaltools.h>
#include <tvb-cpp/tools/eigen/eigen.h>

#include "filter.h"

class BandPassFilter : public Filter {
    double m_TR; //                     # sampling interval
    int k = 2; //                             # 2nd order butterworth filter
    double m_flp = .02; //                        # lowpass frequency of filter
    double m_fhi = 0.1; //                         # highpass

public:
    BandPassFilter(double flp, double fhi, double TR=2.0): m_flp(flp), m_fhi(fhi), m_TR(TR) {}

    tvb::TArray2d apply(const tvb::TArray2d& boldSignal) const override;
};

#endif //TVB_CPP_FILTER_H
