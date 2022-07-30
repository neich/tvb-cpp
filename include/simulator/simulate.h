//
// Created by imartin on 26-Aug-20.
//

#ifndef TVB_CPP_SIMULATE_H
#define TVB_CPP_SIMULATE_H

#include <vector>

#include <simulator/model.h>
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

#include <simulator/integrator.h>
#include <simulator/coupling.h>
#include <datatypes/connectivity.h>
#include <simulator/monitor.h>

namespace tvb {

    enum SIM_MODE {
        SIM_FIXED = 0,
        SIM_ADAPTIVE = 1
    };

    inline
    tvb::TArray2d stateTrackToMatrix(const StateTrack& signal, int axis= 0) {
        assert(signal.m_states.size() > 0);
        TArray2d result(signal.m_states[0].rows(), signal.m_states.size());
        for (unsigned t = 0; t < signal.m_states.size(); ++t)
            result.col(t) = signal.m_states[t].col(axis);
        return result;
    }

    class SimConfig {
    private:
        Model *m_model;
        Integrator *m_integrator;
        Monitor *m_monitor;
        Connectivity *m_connectivity;
        Coupling *m_coupling;
        double m_start_time;
        double m_end_time;
        double m_dt;
        SIM_MODE m_sim_mode;
        double m_delta_integration;
        int m_n_iterations;
        int m_svar_integration;
        int m_sampling_rate;

        int m_n_nodes;
    public:
        SimConfig() :
            m_model(NULL),
            m_integrator(NULL),
            m_connectivity(NULL),
            m_coupling(NULL),
            m_sim_mode(SIM_FIXED),
            m_delta_integration(0.0),
            m_n_iterations(1),
            m_svar_integration(0),
            m_sampling_rate(1)
            {}

        bool is_configured() const {
            return m_model != NULL && m_connectivity != NULL && m_coupling != NULL && m_integrator != NULL;
        }

        void setModel(Model *model) { m_model = model; }

        void setConnectivity(Connectivity *connectivity) {
            m_connectivity = connectivity;
            m_n_nodes = connectivity->weights().rows();
        }

        void setIntegrator(Integrator *integrator) { m_integrator = integrator; }

        void setCoupling(Coupling *coupling) { m_coupling = coupling; }

        void setSamplingRate(int samplingRate) {
            m_sampling_rate = samplingRate;
        }

        Model *model() {
            return m_model;
        }

        Connectivity *connectivity() {
            return m_connectivity;
        }

        Integrator *integrator() {
            return m_integrator;
        }

        Monitor *monitor() {
            return m_monitor;
        }

        Coupling *coupling() {
            return m_coupling;
        }

        int n_nodes() const { return m_n_nodes; }

        double start_time() const {
            return m_start_time;
        }

        double end_time() const {
            return m_end_time;
        }

        double dt() const {
            return m_dt;
        }

        SIM_MODE sim_mode() const { return m_sim_mode; }

        int svar_index() const { return m_svar_integration; }

        double delta_integration() const { return m_delta_integration; }

        int num_iterations() const { return m_n_iterations; }

        int samplingRate() const {
            return m_sampling_rate;
        }

        void setIntegrationInterval(double start_time, double end_time) {
            m_start_time = start_time;
            m_end_time = end_time;
        }

        void setTimeDelta(double dt) {
            m_dt = dt;
        }

        void setMonitor(tvb::Monitor *pMonitor) {
            m_monitor = pMonitor;
        }
    };

    StateTrack simulate(tvb::SimConfig& sim_config);

}


#endif //TVB_CPP_SIMULATE_H
