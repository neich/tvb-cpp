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

#ifndef TVB_CPP_FILTER_DESIGN_H
#define TVB_CPP_FILTER_DESIGN_H

#include <tuple>

#include "tvb-cpp/definitions.h"

typedef std::tuple<tvb::TArray1d, tvb::TArray1dc, tvb::Float> ZeroPoleGain;

ZeroPoleGain iirfilter_zpk(int N,
                           const tvb::TArray1d &Wn,
                           float rp = 0.0,
                           float rs = 0.0,
                           const std::string &btype = "band",
                           bool analog = false,
                           const std::string &ftype = "butter",
                           float fs = 0.0);

std::pair <tvb::TArray1d, tvb::TArray1d> iirfilter_ba(int N,
                                                      const tvb::TArray1d &Wn,
                                                      float rp = 0.0,
                                                      float rs = 0.0,
                                                      const std::string &btype = "bandpass",
                                                      bool analog = false,
                                                      const std::string &ftype = "butter",
                                                      float fs = 0.0);


#endif //TVB_CPP_FILTER_DESIGN_H
