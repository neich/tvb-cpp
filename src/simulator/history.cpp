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


#include <simulator/history.h>

using namespace tvb;

HistoryDense::QResult HistoryDense::query(int step) const {
//            Matrixi time_idx(m_idelays.rows(), m_idelays.cols());
//            tvb::transform(m_idelays, time_idx, (step - 1 - boost::phoenix::placeholders::arg1 + m_ntime) % m_ntime);
    // Matrixd delayed(m_idelays.rows(), m_idelays.cols());
    std::vector<Matrixd> delayed(m_nnodes);
    std::fill_n(delayed.begin(), m_nnodes, tvb::Matrixd(m_nvars, m_nnodes));
    Matrixd current(m_idelays.rows(), m_cvars.size());

    int c_idx = (step - 1) % m_ntime;
    current = m_buffer(c_idx);


    for (unsigned inode = 0; inode < m_idelays.cols(); ++inode)
        for (unsigned inode_from = 0; inode_from < m_idelays.rows(); ++inode_from) {
            int idx = (step - 1 - m_idelays(inode_from, inode) + m_ntime) % m_ntime;
            delayed[inode].col(inode_from) = m_buffer(idx).col(inode_from);
        }

    return History::QResult(current, delayed);
}