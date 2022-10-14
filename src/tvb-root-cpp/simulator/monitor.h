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

#ifndef TVB_CPP_MONITOR_H
#define TVB_CPP_MONITOR_H

#include <tvb-root-cpp/definitions.h>
#include <tvb-root-cpp/simulator/model.h>

namespace tvb {
// Abstract base class for monitor implementations.

    class Monitor {
    public:

        typedef struct {
            float time;
            TArray2d record;
        } Record;

    protected:
        float m_period;
        float m_dt;
        float m_start_time = 0.0;

        int m_istep;
        std::vector<int> m_vars_of_interest;
        std::vector<Record> m_records;
        //        variables_of_interest = NArray(
        //                dtype=int,
        //                label="Model variables to watch",  # order=11,
        //                doc=("Indices of model's variables of interest (VOI) that this monitor should record. "
        //                     "Note that the indices should start at zero, so that if a model offers VOIs V, W and "
        //                     "V+W, and W is selected, and this monitor should record W, then the correct index is 0."),
        //                required=False)

        void init(float period, float dt, const std::vector<int>& voi) {
            m_period = period;
            m_dt = dt;
            m_vars_of_interest = voi;
            m_istep = lround(m_period/m_dt);
        }

    public:
        Monitor(float period, float dt, const std::vector<int>& voi) {
            init(period, dt, voi);
        }

        Monitor(float dt, const std::vector<int>& voi): m_dt(dt), m_vars_of_interest(voi) {
            init(1.0, dt, voi);
        }

        virtual ~Monitor() = default;

        void setStartTime(float st) { m_start_time = st; }

        virtual void from_records(const std::vector<Record>& from, std::vector<Record>& to) {
            throw std::runtime_error("Method from_records() not implemented");
        }

        void record(int step, const State &observed) {
            this->sample(step, observed);
        }

        virtual void sample(int step, const State &state) = 0;

        const std::vector<Record>& getRecords() const {
            return m_records;
        }

        [[nodiscard]] TArray2d voi2Array(int index) const {
            int N = getRecords()[0].record.rows();
            TArray2d result(N, getRecords().size());
            for (int i = 0; i < getRecords().size(); ++i)
                result.col(i) = getRecords()[i].record.col(index);
            return result;
        }

    };

    class Raw : public Monitor {

    public:
        Raw(float dt, const std::vector<int>& voi): Monitor(dt, voi) {}

        void sample(int step, const State &state) override {
            m_records.push_back(Record{m_start_time + step * m_dt, state(Eigen::all, m_vars_of_interest)});
        }
    };

    class RawSubSample : public Monitor {
        int m_every_n;

    public:
        explicit RawSubSample(float period, float dt, const std::vector<int>& voi): Monitor(period, dt, voi) {}

        void sample(int step, const State &state) override {
            if (step % m_istep == 0)
                m_records.push_back(Record{m_start_time + step * m_dt, state(Eigen::all, m_vars_of_interest)});
        }

    };

    class AverageSubSample : public Monitor {
        std::vector<State> m_buffer;

    public:
        explicit AverageSubSample(int n_nodes, float period, float dt, const std::vector<int>& voi): Monitor(period, dt, voi) {
            m_buffer.resize(m_istep);
            std::fill_n(m_buffer.begin(), m_istep, TArray2d::Zero(n_nodes, voi.size()));
        }

        void sample(int step, const State &state) override {
            m_buffer[step % m_istep] = state(Eigen::all, m_vars_of_interest);
            if (step % m_istep == 0) {
                TArray2d result = TArray2d::Zero(state.rows(), m_vars_of_interest.size());
                for (auto &s: m_buffer)
                    result += s;
                result /= m_istep;
                m_records.push_back(Record{m_start_time + step * m_dt, result});
            }
        }

    };
}

#endif //TVB_CPP_MONITOR_H
