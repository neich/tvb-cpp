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

#ifndef TVB_CPP_CONNECTIVITY_H
#define TVB_CPP_CONNECTIVITY_H

#include <cassert>
#include <algorithm>

#include <tvb-root-cpp/definitions.h>

namespace tvb {

    class Connectivity {
        int m_num_nodes;
        TArray2d m_weights;
        TArray2d m_tract_lengths;
        double m_speed;

        TArray2d m_delays;

    public:

        Connectivity(int n_nodes = 0) : m_num_nodes(n_nodes) {}

        Connectivity(const TArray2d &weights, const TArray2d &tract_lengths, double speed) {
            assert(weights.rows() == weights.cols());
            assert(tract_lengths.rows() == tract_lengths.cols());
            assert(tract_lengths.rows() == weights.rows());
            m_num_nodes = weights.rows();
            this->setWeights(weights);
            this->setTractLengths(tract_lengths);
            this->setSpeed(speed);
            m_delays = m_tract_lengths / m_speed;

        }

        void setWeights(const TArray2d &weights) {
            assert(weights.rows() == m_num_nodes);
            assert(weights.cols() == m_num_nodes);
            m_weights = weights;
        }

        void setTractLengths(const TArray2d &tract_lengths) {
            assert(tract_lengths.rows() == m_num_nodes);
            assert(tract_lengths.cols() == m_num_nodes);
            m_tract_lengths = tract_lengths;
        }

        const TArray2d &weights() const { return m_weights; }

        const TArray2d &delays() const { return m_delays; }

        void setSpeed(double speed) { m_speed = speed; }
    };

}

#endif //TVB_CPP_CONNECTIVITY_H
