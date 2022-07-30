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
    protected:
        int m_n_nodes;
        int m_n_voi;

        double m_period;
        //        period = Float(
        //                label="Sampling period (ms)",  # order = 10
        //        default=0.9765625,  # ms. 0.9765625 => 1024Hz #ms, 0.5 => 2000Hz
        //        doc="""Sampling period in milliseconds, must be an integral multiple
        //        of integration-step size. As a guide: 2048 Hz => 0.48828125 ms ;
        //        1024 Hz => 0.9765625 ms ; 512 Hz => 1.953125 ms.""")

        std::vector<int> m_vars_of_interest;
        //        variables_of_interest = NArray(
        //                dtype=int,
        //                label="Model variables to watch",  # order=11,
        //                doc=("Indices of model's variables of interest (VOI) that this monitor should record. "
        //                     "Note that the indices should start at zero, so that if a model offers VOIs V, W and "
        //                     "V+W, and W is selected, and this monitor should record W, then the correct index is 0."),
        //                required=False)

        int m_istep;
        double m_dt;
        //        voi = None
        //        _stock = numpy.empty([])

    public:
        Monitor() : m_n_voi(1), m_period(0.9765625), m_vars_of_interest({0})  {}


        double period() const { return m_period; }

        //def __str__(self):
        //clsname = self.__class__.__name__
        //return '%s(period=%f, voi=%s)' % (clsname, self.period, self.variables_of_interest.tolist())

//        virtual void config_for_sim(const SimConfig &sim_config, const std::vector<int> &voi) {
//            m_n_nodes = sim_config.n_nodes();
//            m_n_voi = voi.size();
//            //"""Configure monitor for given simulator.
//            //
//            //Grab the Simulator's integration step size. Set the monitor's variables
//            //        of interest based on the Monitor's 'variables_of_interest' attribute, if
//            //it was specified, otherwise use the 'variables_of_interest' specified
//            //for the Model. Calculate the number of integration steps (isteps)
//            //between returns by the record method. This method is called from within
//            //the the Simulator's configure() method.
//            //
//            //"""
//            m_dt = sim_config.dt();
//            m_istep = int(round(m_period / m_dt));
//            m_vars_of_interest = voi;
//        }
//
//        virtual void config_for_sim(int N, double dt, const std::vector<int> &voi) {
//            m_n_nodes = N;
//            m_n_voi = voi.size();
//            //"""Configure monitor for given simulator.
//            //
//            //Grab the Simulator's integration step size. Set the monitor's variables
//            //        of interest based on the Monitor's 'variables_of_interest' attribute, if
//            //it was specified, otherwise use the 'variables_of_interest' specified
//            //for the Model. Calculate the number of integration steps (isteps)
//            //between returns by the record method. This method is called from within
//            //the the Simulator's configure() method.
//            //
//            //"""
//            m_dt = dt;
//            m_istep = int(round(m_period / m_dt));
//            m_vars_of_interest = voi;
//        }


        void record(int step, const State &observed) {
            //"""Record a sample of the observed state at given step.
            //
            //This is a final method called by the simulator to obtain samples from a
            //monitor instance. Monitor subclasses should not override this method, but
            //        rather implement the `sample` method.
            //
            //"""
            this->sample(step, observed);
        }

        virtual void sample(int step, const State &state) = 0;
        //"""
        //This method provides monitor output, and should be overridden by subclasses.
        //
        //"""


        //    void create_time_series(const Connectivity* connectivity) {
        //"""
        //Create a time series instance that will be populated by this monitor
        //:param surface: if present a TimeSeriesSurface is returned
        //:param connectivity: if present a TimeSeriesRegion is returned
        //        Otherwise a plain TimeSeries will be returned
        //"""
        //        if surface
        //        is
        //        not None:
        //        return TimeSeriesSurface(surface = surface.region_mapping_data.surface,
        //                                 sample_period = self.period,
        //                                 title = 'Surface ' + self.__class__.__name__)
        //        if connectivity
        //        is
        //        not None:
        //        return TimeSeriesRegion(connectivity = connectivity,
        //                                region_mapping = region_map,
        //                                region_mapping_volume = region_volume_map,
        //                                sample_period = self.period,
        //                                title = 'Regions ' + self.__class__.__name__)
        //
        //        return TimeSeries(sample_period = self.period,
        //                          title = ' ' + self.__class__.__name__)
        //    }
//        virtual StateTrack apply(const std::vector<double>& times, const std::vector<State>& states) {
//            StateTrack result;
//            for (unsigned step = 0; step < times.size(); ++step) {
//                auto [time, sample] = this->sample(step+1, states[step]);
//                if (time >= 0.0) result.push(sample, time);
//            }
//            return result;
//        }
//
//        virtual TArray2d apply(const TArray2d& signal) {
//            TArray2d result(signal.rows(), signal.cols());
//            for (unsigned step = 0; step < signal.size(); ++step) {
//                auto [time, sample] = this->sample(step+1, signal.col(step));
//                result.col(step) = sample;
//            }
//            return result;
//        }

    };

    class Raw : public Monitor {
        std::vector<State> m_records;

    public:
        Raw() = default;

        void sample(int step, const State &state) override {
            m_records.push_back(state);
        }
    };

    class RawSubSample : public Monitor {
        std::vector<State> m_records;
        int m_every_n;

    public:
        explicit RawSubSample(int every_n): m_every_n(every_n) {}

        void sample(int step, const State &state) override {
            if (step % m_every_n == 0)
                m_records.push_back(state);
        }

        [[nodiscard]] const std::vector<State>& getRecords() const { return m_records; }
    };


}

#endif //TVB_CPP_MONITOR_H
