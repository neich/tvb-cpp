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

#ifndef TVB_CPP_HISTORY_H
#define TVB_CPP_HISTORY_H

#include <definitions.h>
#include <simulator/model.h>

namespace tvb {

    class History {
    public:
        typedef typename std::pair<Matrixd, std::vector<Matrixd>> QResult;

    protected:
        Matrixd m_weights;
        Matrixd m_delays;
        std::vector<int> m_cvars;
        int m_nvars;
        int m_nnodes;

    public:

        History() {}

        History(int n_nodes, int n_cvar) : m_weights(n_nodes, n_nodes), m_cvars(n_cvar) {}

        History(const Matrixd &weights, const Matrixd &delays, const std::vector<int> &cvars) {
            m_weights = weights;
            m_delays = delays;
            m_cvars = cvars;
        }

        int num_nodes() const { return m_nnodes; }

        int num_vars() const { return m_nvars; }

        const std::vector<int>& c_vars() const { return m_cvars; }

        [[nodiscard]] const Matrixd &weights() const { return m_weights; }

        Matrixd &delays() { return m_delays; }

        std::vector<int> &cvars() { return m_cvars; }

        virtual void init(double dt, const State& init_state) {
            m_nvars = init_state.cols();
            m_nnodes = init_state.rows();
        }

        [[nodiscard]] virtual QResult query(int step) const = 0;

        virtual void update(int step, int n_reg, const State& state) = 0;

    };

    class HistoryDense : public History {
    public:
        typedef typename Eigen::Matrix<Matrixd, Eigen::Dynamic, 1> Buffer; // (time, vars, nodes)

    private:
        Matrixi m_idelays;
        int m_ntime = 1;

        Buffer m_buffer;

    public:

        HistoryDense(): History() {}

        // HistoryDense(int n_nodes, int n_cvar) : History(n_nodes, n_cvar) {}

        HistoryDense(const Matrixd &weights, const Matrixd &delays, const std::vector<int> &cvars): History(weights, delays, cvars) {}

        void init(double dt, const State& init_state) override {
            History::init(dt, init_state);
            m_idelays = Matrixi(m_weights.rows(), m_weights.cols());
            // tvb::transform(m_delays, m_idelays, rint(boost::phoenix::placeholders::arg1 / dt));
            tvb::transform(m_delays, m_idelays, [&dt](double arg1) { return int(lround(arg1 / dt)); });
            m_ntime = m_idelays.maxCoeff() + 1;
            m_buffer = Buffer(m_ntime);
            tvb::fill(m_buffer, init_state);
        }

        [[nodiscard]] QResult query(int step) const override;

        void update(int step, int n_reg, const State& state) override {
            m_buffer(step % m_ntime) = state(Eigen::all, m_cvars);
        }

    };

    class HistoryNoDelays : public History {
    private:

        tvb::Matrixd m_buffer;

    public:

        HistoryNoDelays(): History() {}

        HistoryNoDelays(int n_nodes, int n_cvar) : History(n_nodes, n_cvar) {}

        HistoryNoDelays(const Matrixd &weights, const Matrixd &delays, const std::vector<int> &cvars): History(weights, delays, cvars) {}

        void init(double dt, const State& init_state) override {
            History::init(dt, init_state);
            m_buffer = init_state(Eigen::all, m_cvars);
        }

        QResult query(int step) const override {
            std::vector<tvb::Matrixd> delayed(m_nnodes);
            std::fill_n(delayed.begin(), m_nnodes, m_buffer.transpose());
            return History::QResult(m_buffer, delayed);
        }

        void update(int step, int n_reg, const State& state) override {
            m_buffer = state(Eigen::all, m_cvars);
        }

    };


}

#endif //TVB_CPP_HISTORY_CPP
