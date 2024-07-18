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

#ifndef TVB_CPP_COUPLING_H
#define TVB_CPP_COUPLING_H

#include <tvb-cpp/definitions.h>
#include "history.h"

namespace tvb {

    class Coupling {
    protected:
        TArray2d m_weights;
        TArray2d m_delays;
        std::vector<int> m_cvars;
        int m_nvars;
        int m_nnodes;

    public:

        Coupling(const TArray2d &weights, const TArray2d &delays, const std::vector<int> &cvars) {
            m_weights = weights;
            m_delays = delays;
            m_cvars = cvars;
            m_nnodes = weights.cols();
            m_nvars = cvars.size();
        }

        virtual void init(double dt, const State &init_state) = 0;

        virtual TArray2d couple(int step) const = 0;

        virtual std::vector<tvb::TArray2d> pre(const TArray2d &x_i, std::vector<tvb::TArray2d> &x_j) const {
            return x_j;
        }

        virtual TArray2d post(const TArray2d &gx) const {
            return gx;
        }

        virtual void update(int step, const State &state) = 0;

    };

    class CouplingNoDelays : public Coupling {
    protected:
        int m_ntime;
        double m_dt;
        double m_scale;

        TArray2d m_buffer;

    public:
        CouplingNoDelays(const TArray2d &weights, const TArray2d &delays, const std::vector<int> &cvars) : Coupling(
                weights, delays, cvars), m_scale(1.0) {};

        void init(double dt, const State &init_state) override {
            m_dt = dt;
            m_buffer = TArray2d(this->m_nnodes, this->m_nvars);
            // tvb::transform(m_delays, m_idelays, rint(boost::phoenix::placeholders::arg1 / dt));
            this->update(0, init_state);
        }

        virtual TArray2d couple(int step) const override {
            TArray2d result(m_nnodes, m_nvars);

            for (int c = 0; c < m_cvars.size(); ++c) {
                result.col(c) = m_weights.matrix() * this->m_buffer.col(c).matrix();
            }

            return result;
        }

        virtual void update(int step, const State &state) override {
            for (int c = 0; c < m_nvars; ++c)
                m_buffer.col(c) = state.col(m_cvars[c]);
        }

    };


    class CouplingLinearDense : public Coupling {
    protected:
        TArray2di m_idelays;
        int m_ntime;
        double m_dt;

        std::vector<TArray2d> m_buffer;

    public:
        CouplingLinearDense(const TArray2d &weights, const TArray2d &delays, const std::vector<int> &cvars) : Coupling(
                weights, delays, cvars) {};

        void init(double dt, const State &init_state) override {
            m_dt = dt;
            m_idelays = TArray2di(m_weights.rows(), m_weights.cols());
            // tvb::transform(m_delays, m_idelays, rint(boost::phoenix::placeholders::arg1 / dt));
            tvb::transform(m_delays, m_idelays, [&dt](double arg1) { return int(lround(arg1 / dt)); });
            m_ntime = m_idelays.maxCoeff() + 1;
            m_buffer = std::vector<TArray2d>(m_cvars.size());
            std::fill_n(m_buffer.begin(), m_cvars.size(), TArray2d::Zero(this->m_nnodes, m_ntime));
            this->update(0, init_state);
        }


        virtual TArray2d couple(int step) const override;

        virtual void update(int step, const State &state) override {
            for (int c = 0; c < m_cvars.size(); ++c)
                m_buffer[c].col(step % m_ntime) = state.col(c);
        }

    };

    class CouplingLinearSparse : public Coupling {
    protected:
        TArray2di m_idelays;
        int m_ntime = 1;
        double m_dt;
        double m_scale;

        std::vector<TArray2d> m_buffer;

        TArray1di m_index_sizes;
        std::vector<TArray1d> m_wsparse;
        std::vector<TArray1di> m_dsparse;
        std::vector<std::vector<std::vector<Float*>>> m_pbuffer;


    public:
        CouplingLinearSparse(const TArray2d &weights, const TArray2d &delays, const std::vector<int> &cvars) : Coupling(
                weights, delays, cvars), m_scale(1.0) {};

        void init(double dt, const State &init_state) override;


        TArray2d couple(int step) const override;

        void update(int step, const State &state) override {
            for (int c = 0; c < m_cvars.size(); ++c)
                m_buffer[c].col(step % m_ntime) = state.col(m_cvars[c]);
        }

        void setScale(double scale) {
            m_scale = scale;
        }
    };

}

#endif //TVB_CPP_COUPLING_CPP
