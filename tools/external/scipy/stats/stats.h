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

#ifndef TVB_CPP_MATH_H
#define TVB_CPP_MATH_H


#include <definitions.h>

using namespace tvb;

double corrcoef(const TArray1d &x, const TArray1d &y);

std::pair<double, double> ks_2samp(const TArray1d &data1,
                                   const TArray1d &y,
                                   const std::string& alternative="two-sided",
                                   std::string mode="auto");

#endif //TVB_CPP_MATH_H
