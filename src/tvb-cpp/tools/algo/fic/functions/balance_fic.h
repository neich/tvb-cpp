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

#ifndef TVB_CPP_BALANCE_FIC_H
#define TVB_CPP_BALANCE_FIC_H

#include <tvb-cpp/definitions.h>
#include <tvb-cpp/simulator/model.h>
#include <tvb-cpp/simulator/simulate.h>


struct OptResult {
    std::vector<tvb::State> m_states;
    std::vector<double> m_times;
    tvb::TArray1d m_Jis;
};


std::tuple<int, double, tvb::TArray1d>  optimize_fic(tvb::SimConfig& sim_config, int voi, float value_base);
std::tuple<bool, float, float, double, tvb::TArray1d>  optimize_fic_Herzog(tvb::SimConfig& sim_config, int voi, float value_base, float a = 1.0, float b = 0.5);

#endif //TVB_CPP_BALANCE_FIC_H
