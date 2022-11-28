# -*- coding: utf-8 -*-
#
#
# TheVirtualBrain-Scientific Package. This package holds all simulators, and
# analysers necessary to run brain-simulations. You can use it stand alone or
# in conjunction with TheVirtualBrain-Framework Package. See content of the
# documentation-folder for more details. See also http://www.thevirtualbrain.org
#
# (c) 2012-2022, Baycrest Centre for Geriatric Care ("Baycrest") and others
#
# This program is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this
# program.  If not, see <http://www.gnu.org/licenses/>.
#
#
#   CITATION:
# When using The Virtual Brain for scientific publications, please cite it as follows:
#
#   Paula Sanz Leon, Stuart A. Knock, M. Marmaduke Woodman, Lia Domide,
#   Jochen Mersmann, Anthony R. McIntosh, Viktor Jirsa (2013)
#       The Virtual Brain: a simulator of primate brain network dynamics.
#   Frontiers in Neuroinformatics (7:10. doi: 10.3389/fninf.2013.00010)
#
#

"""
C++ backend which uses templating to generate simulation code.

.. moduleauthor:: Ignacio Martín <ignacio.martin@udg.edu>

"""

import os
import importlib
import numpy as np

from .templates import MakoUtilMix

import autopep8
from tvb.simulator.lab import *
from tvb.basic.neotraits.api import NArray

import tvbcpp


class CppBackend(MakoUtilMix):

    def build_py_func(self, template_source, content, name='kernel', print_source=False,
            modname=None):
        "Build and retrieve one or more Python functions from template."
        source = self.render_template(template_source, content)
        source = autopep8.fix_code(source)
        if print_source:
            print(self.insert_line_numbers(source))
        if modname is not None:
            return self.eval_module(source, name, modname)
        else:
            return self.eval_source(source, name)

    def eval_source(self, source, name, print_source=False):
        globals_ = {}
        try:
            exec(source, globals_)
        except Exception as exc:
            if not print_source:
                print(self._insert_line_numbers(source))
            raise exc
        fns = [globals_[n] for n in name.split(',')]
        return fns[0] if len(fns)==1 else fns

    def eval_module(self, source, name, modname):
        here = os.path.abspath(os.path.dirname(__file__))
        genp = os.path.join(here, 'templates', 'generated')
        with open(f'{genp}/{modname}.py', 'w') as fd:
            fd.write(source)
        fullmodname = f'tvb.simulator.backend.templates.generated.{modname}'
        mod = importlib.import_module(fullmodname)
        fns = [getattr(mod,n) for n in name.split(',')]
        return fns[0] if len(fns)==1 else fns

    def check_compatibility(self, sim):
        def check_choices(val, choices):
            if not isinstance(val, choices):
                raise NotImplementedError(
                    "Unsupported simulator component. Given: {}\nExpected one of: {}".format(val, choices))

        # monitors
        assert len(sim.monitors) == 1
        for m in sim.monitors:
            check_choices(m, (monitors.Raw, monitors.TemporalAverage))  # looks like this is the only one in the python API now

        # integrators
        check_choices(
            sim.integrator,
            (  # uncomment when tested
                integrators.EulerStochastic,
                integrators.EulerDeterministic,
                # integrators.HeunStochastic,
                # integrators.HeunDeterministic,
            )
        )
        # models
        check_choices(
            sim.model,
            (  # uncomment when tested
                # models.MontbrioPazoRoxin,
                models.ReducedWongWangExcInh,
                models.ZerlautAdaptationFirstOrder,
                models.ZerlautAdaptationSecondOrder
            )
        )

        # coupling
        check_choices(sim.coupling, coupling.SparseCoupling)
        # surface
        if sim.surface is not None:
            raise NotImplementedError("Surface simulation not supported.")
        # stimulus
        if sim.stimulus is not None:
            raise NotImplementedError("Stimulation not supported.")

    def _build_cpp_sim(self, sim):
        tvbcpp.set_weights(
            sim.connectivity.weights
        )
        tvbcpp.set_lenghts(
            sim.connectivity.tract_lengths,
            sim.connectivity.speed
        )

        print(sim.integrator.dt)

        tvbcpp.set_integrator_es(
            sim.integrator.dt,
            np.append(sim.integrator.noise.nsig.squeeze(), [0., 0.])
            # !!! this is wrong in many cases
        )

        # TODO case map between the simulators
        tvbcpp.set_model(type(sim.model).__name__)
        for param_name in type(sim.model).declarative_attrs:
            param_def = getattr(type(sim.model), param_name)
            if not isinstance(param_def, NArray) or not param_def.dtype == np.float:
                continue
            param_value = getattr(sim.model, param_name)
            tvbcpp.set_model_parameter(param_name, param_value)

        for m in sim.monitors:
            if isinstance(m, monitors.Raw):
                tvbcpp.add_raw_monitor(m.period, m.voi)
            elif isinstance(m, monitors.TemporalAverage):
                tvbcpp.add_average_monitor(m.period, m.voi)

    def run_sim(self, sim, nstep=None, simulation_length=None, chunksize=100000, compatibility_mode=False):
        assert nstep is not None or simulation_length is not None or sim.simulation_length is not None

        if simulation_length is None:
            simulation_length = sim.simulation_length

        self.check_compatibility(sim)
        self._build_cpp_sim(sim)
        d_raw = tvbcpp.run_sim(0, simulation_length)

        nstep = int(np.ceil(simulation_length / sim.integrator.dt))
        t_raw = np.arange(nstep) * sim.integrator.dt
        return (t_raw, d_raw[:, :, :, np.newaxis]),