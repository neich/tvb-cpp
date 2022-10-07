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

#include <simulator/coupling.h>

using namespace tvb;

TArray2d CouplingLinearDense::couple(int step) const {
    TArray2d result(m_nnodes, m_nvars);

    for (unsigned cv = 0; cv < m_nvars; ++cv) {
        TArray2d delayed(m_nnodes, m_nnodes);

        for (unsigned inode = 0; inode < m_nnodes; ++inode) {
            for (unsigned inode_from = 0; inode_from < m_nnodes; ++inode_from) {
                int idx = (step - 1 - m_idelays(inode, inode_from) + m_ntime) % m_ntime;
                delayed(inode_from, inode) = m_buffer[cv](inode_from, idx);
            }
        }

        for (unsigned inode = 0; inode < m_nnodes; ++inode) {
            result(inode, cv) = (delayed.col(inode) * m_weights.col(inode)).sum();
        }
    }

    return result;
}

void CouplingLinearSparse::init(double dt, const State &init_state) {
    m_dt = dt;
    m_idelays = TArray2di(m_weights.rows(), m_weights.cols());
    // tvb::transform(m_delays, m_idelays, rint(boost::phoenix::placeholders::arg1 / dt));
    tvb::transform(m_delays, m_idelays, [&dt](Float arg1) { return int(lround(arg1 / dt)); });
    m_ntime = m_idelays.maxCoeff() + 1;
    m_buffer = std::vector<TArray2d>(m_cvars.size());
    std::fill_n(m_buffer.begin(), m_cvars.size(), TArray2d::Zero(m_nnodes, m_ntime));
    this->update(0, init_state);

    m_index_sizes.resize(m_nnodes);
    m_wsparse.resize(m_nnodes);
    m_dsparse.resize(m_nnodes);

    std::vector<int> indices;
    for (int n = 0; n < m_nnodes; ++n) {
        int nsize = 0;
        for (int nn = 0; nn < m_nnodes; ++nn) {
            if (m_weights(n, nn) > 0.0) {
                indices.push_back(nn);
                nsize++;
            }
        }
        m_index_sizes[n] = nsize;
    }

    int i = 0;
    for (int n = 0; n < m_nnodes; ++n) {
        int nsize = m_index_sizes[n];
        m_wsparse[n] = TArray1d(nsize);
        m_dsparse[n] = TArray1di(nsize);
        for (int in = 0; in < nsize; ++in) {
            m_wsparse[n][in] = m_weights(n, indices[i]);
            m_dsparse[n][in] = m_idelays(n, indices[i]);
            i++;
        }
    }

    m_pbuffer.clear();
    for (int cv = 0; cv < m_nvars; ++cv) {
        m_pbuffer.emplace_back();
        int index = 0;
        for (int n = 0; n < m_nnodes; ++n) {
            m_pbuffer.back().emplace_back();
            std::vector<Float *> &md = m_pbuffer.back().back();
            int nsize = m_index_sizes[n];
            const TArray1di &dsparse = m_dsparse[n];
            for (int step = 0; step < m_ntime; ++step) {
                for (int nn = 0; nn < nsize; ++nn) {
                    int idx = (step - 1 - dsparse[nn] + m_ntime) % m_ntime;
                    Float *fp = &m_buffer[cv](indices[index+nn], idx);;
                    md.push_back(fp);
                }
            }
            index += nsize;
        }
    }

}

TArray2d CouplingLinearSparse::couple(int step) const {
    TArray2d result(m_nnodes, m_nvars);

    int step_cycle = step % m_ntime;


    TArray1d delayed(m_nnodes);
    for (unsigned c = 0; c < m_nvars; ++c) {
//        if (tvb::isnan(m_buffer[c]))
//            std::cout << "Ein"; // TODO debug!
        for (unsigned n = 0; n < m_nnodes; ++n) {
            int nsize = m_index_sizes[n];
            const TArray1d& wsparse = m_wsparse[n];
            const std::vector<Float*>& bpointers = m_pbuffer[c][n];
            Float sum = 0.0;
            for (int s = 0; s < nsize; ++s) {
                Float* fp = bpointers[step_cycle * nsize + s];
                sum += *fp * wsparse[s];
//                if (std::isnan(sum))
//                    std::cout << "ein";
            }
            result(n, c) = sum;
        }
    }

    return result;

}



