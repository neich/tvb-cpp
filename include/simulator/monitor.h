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

#include <definitions.h>
#include <simulator/model.h>

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

        void record(int step, const State &observed) {
            this->sample(step, observed);
        }

        virtual void sample(int step, const State &state) = 0;

        const std::vector<Record>& getRecords() const {
            return m_records;
        }

    };

    class Raw : public Monitor {

    public:
        Raw(float dt, const std::vector<int>& voi): Monitor(dt, voi) {}

        void sample(int step, const State &state) override {
            m_records.push_back(Record{step * m_dt, state(Eigen::all, m_vars_of_interest)});
        }
    };

    class RawSubSample : public Monitor {
        int m_every_n;

    public:
        explicit RawSubSample(float period, float dt, const std::vector<int>& voi): Monitor(period, dt, voi) {}

        void sample(int step, const State &state) override {
            if (step % m_istep == 0)
                m_records.push_back(Record{step * m_dt, state(Eigen::all, m_vars_of_interest)});
        }

    };
}

#endif //TVB_CPP_MONITOR_H
