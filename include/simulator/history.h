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
        typedef typename std::pair<TArray2d, std::vector<TArray2d>> QResult;

    protected:
        TArray2d m_weights;
        TArray2d m_delays;
        std::vector<int> m_cvars;
        int m_nvars;
        int m_nnodes;

    public:

        History() {}

        History(int n_nodes, int n_cvar) : m_weights(n_nodes, n_nodes), m_cvars(n_cvar) {}

        History(const TArray2d &weights, const TArray2d &delays, const std::vector<int> &cvars) {
            m_weights = weights;
            m_delays = delays;
            m_cvars = cvars;
        }

        int num_nodes() const { return m_nnodes; }

        int num_vars() const { return m_nvars; }

        const std::vector<int>& c_vars() const { return m_cvars; }

        [[nodiscard]] const TArray2d &weights() const { return m_weights; }

        TArray2d &delays() { return m_delays; }

        std::vector<int> &cvars() { return m_cvars; }

        virtual void init(double dt, const State& init_state) {
            m_nvars = init_state.cols();
            m_nnodes = init_state.rows();
        }

        [[nodiscard]] virtual QResult query(int step) const = 0;

        virtual void update(int step, const State& state) = 0;

        virtual const TArray2d &getCurrent() const = 0;
    };

    class HistoryDense : public History {
    public:

    private:
        TArray2di m_idelays;
        int m_ntime = 1;

        std::vector<TArray2d> m_buffer;

    public:

        HistoryDense(): History() {}

        // HistoryDense(int n_nodes, int n_cvar) : History(n_nodes, n_cvar) {}

        HistoryDense(const TArray2d &weights, const TArray2d &delays, const std::vector<int> &cvars): History(weights, delays, cvars) {}

        void init(double dt, const State& init_state) override {
            History::init(dt, init_state);
            m_idelays = TArray2di(m_weights.rows(), m_weights.cols());
            // tvb::transform(m_delays, m_idelays, rint(boost::phoenix::placeholders::arg1 / dt));
            tvb::transform(m_delays, m_idelays, [&dt](double arg1) { return int(lround(arg1 / dt)); });
            m_ntime = m_idelays.maxCoeff() + 1;
            m_buffer = std::vector<TArray2d>(m_cvars.size());
//            std::fill_n(m_buffer, m_cvars.size(), TArray2d(this->num_nodes(), m_ntime));
            this->update(0, init_state);
        }

        [[nodiscard]] QResult query(int step) const override;

        void update(int step, const State& state) override {
            for (int c = 0; c < m_cvars.size(); ++c)
                m_buffer[c].col(step % m_ntime) = state.col(c);
        }

        virtual const TArray2d &getCurrent() const override {
            throw "Cannot use getCurrent() with HistoryDense";
        }

    };

    class HistoryNoDelays : public History {
    private:

        tvb::TArray2d m_buffer;

    public:

        HistoryNoDelays(): History() {}

        HistoryNoDelays(const TArray2d &weights, const TArray2d &delays, const std::vector<int> &cvars): History(weights, delays, cvars) {}

        void init(double dt, const State& init_state) override {
            History::init(dt, init_state);
            m_buffer = init_state(Eigen::all, m_cvars);
        }

        [[nodiscard]] QResult query(int step) const override {
            std::vector<tvb::TArray2d> delayed(m_nnodes);
            std::fill_n(delayed.begin(), m_nnodes, m_buffer);
            return {m_buffer, delayed};
        }

        void update(int step, const State& state) override {
            m_buffer = state(Eigen::all, m_cvars);
        }

        [[nodiscard]] const TArray2d &getCurrent() const override {
            return m_buffer;
        }
    };


}

#endif //TVB_CPP_HISTORY_CPP
