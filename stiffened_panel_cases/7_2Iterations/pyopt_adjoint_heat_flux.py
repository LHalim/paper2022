#!/usr/bin/env python
"""
This file is part of the package FUNtoFEM for coupled aeroelastic simulation
and design optimization.

Copyright (C) 2015 Georgia Tech Research Corporation.
Additional copyright (C) 2015 Kevin Jacobson, Jan Kiviaho and Graeme Kennedy.
All rights reserved.

FUNtoFEM is licensed under the Apache License, Version 2.0 (the "License");
you may not use this software except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

from __future__ import print_function

#from os import environ
#environ['CMPLX_MODE'] = "1"
from pyfuntofem.model  import *
from pyfuntofem.driver import *
from pyfuntofem.fun3d_interface import *

from tacs_model import wedgeTACS
from pyoptsparse import SLSQP, Optimization
from mpi4py import MPI
import os
import numpy as np

class wedge_adjoint(object):
    """
    -------------------------------------------------------------------------------
    TOGW minimization
    -------------------------------------------------------------------------------
    """

    def __init__(self, analysis_type='aeroelastic'):
        print('start')

        self.analysis_type = analysis_type

        # cruise conditions
        self.v_inf = 171.5                  # freestream velocity [m/s]
        self.rho = 0.01841                  # freestream density [kg/m^3]
        self.cruise_q = 12092.5527126               # dynamic pressure [N/m^2]
        self.grav = 9.81                            # gravity acc. [m/s^2]
        self.thermal_scale = 0.5 * self.rho * (self.v_inf)**3

        self.maximum_mass = 40.0 
        self.num_tacs_dvs = 112

        # Set up the communicators
        n_tacs_procs = 1

        comm = MPI.COMM_WORLD
        self.comm = comm
        print('set comm')

        world_rank = comm.Get_rank()
        if world_rank < n_tacs_procs:
            color = 55
            key = world_rank
        else:
            color = MPI.UNDEFINED
            key = world_rank
        self.tacs_comm = comm.Split(color,key)
        print('comm misc')
        # Set up the FUNtoFEM model for the TOGW problem
        self._build_model()
        print('built model')
        self.ndv = len(self.model.get_variables())
        print("ndvs: ",self.ndv)
        # instantiate TACS on the master
        solvers = {}
        solvers['flow'] = Fun3dInterface(self.comm,self.model,flow_dt=1.0)#flow_dt=0.1
        solvers['structural'] = wedgeTACS(self.comm,self.tacs_comm,self.model,n_tacs_procs)

        # L&D transfer options
        transfer_options = {'analysis_type': self.analysis_type,
                            'scheme': 'meld', 'thermal_scheme': 'meld'}

        # instantiate the driver
        self.driver = FUNtoFEMnlbgs(solvers,self.comm,self.tacs_comm,0,self.comm,0,transfer_options,model=self.model)

        # Set up some variables and constants related to the problem
        self.cruise_lift   = None
        self.cruise_drag   = None
        self.num_con = 1
        self.mass = None

        self.functions = None
        self.grads = None

        self.var_scale        = np.ones(self.ndv,dtype=TransferScheme.dtype)
        self.struct_tacs = solvers['structural'].assembler

    def _build_model(self):

        thickness = 0.015
        # Build the model
        model = FUNtoFEMmodel('wedge')
        plate = Body('plate', analysis_type=self.analysis_type, group=0,boundary=1)

        for i in range(self.num_tacs_dvs):
            plate.add_variable('structural',Variable('thickness '+ str(i),value=thickness,lower = 0.001, upper = 0.1))

        #plate.add_variable('structural',Variable('thickness',value=thickness,lower = 0.01, upper = 0.1))
        model.add_body(plate)

        steady = Scenario('steady', group=0, steps=2)
        #steady.set_variable('aerodynamic',name='AOA',value=0.0,lower=-15.0,upper=15.0)
        temp = Function('ksfailure',analysis_type='structural') #temperature
        steady.add_function(temp)

        mass = Function('mass',analysis_type='structural',adjoint=False)
        steady.add_function(mass)

        #lift = Function('cl',analysis_type='aerodynamic')
        #steady.add_function(lift)

        #drag = Function('cd',analysis_type='aerodynamic')
        #steady.add_function(drag)

        model.add_scenario(steady)

        self.model = model

    def total_derivative(self):
        self.driver.solve_forward()
        functions = self.model.get_functions()
        
        for index, func in enumerate(functions):
            print('Function %d'%(index), func.value)
        
        self.driver.solve_adjoint()

        grads = self.model.get_function_gradients()
        
        variables = self.model.get_variables()

        for i, func in enumerate(functions):
            for j, var in enumerate(variables):
                print("Grad ", func.name, "Var: ", var.name, " ", grads[i][j])
        
        self.functions = functions
        self.grads = grads
        
        return

    def objFunc(self, xdict):

        tInput = xdict["xvars"]
        self.driver.solve_forward()
        functions = self.model.get_functions()

        funcs = {}
        funcs["obj"] = functions[0].value
        funcs["con"] = functions[1].value

        fail = False

        if self.comm.rank == 0:
            print(tInput) 

        return funcs, fail

    def objGrad(self, xdict, funcs): 

        tInput = xdict["xvars"]

        self.driver.solve_adjoint()
        grads = self.model.get_function_gradients()

        grad1 = np.array(grads[0][:])
        grad2 = np.array(grads[1][:])

        sens = {}
        sens = {
            "obj": {
                "xvars": [grad1]
            },
            "con": {
                "xvars": [grad2]
            },
        }

        fail = False

        return sens, fail

################################################################################
dp = wedge_adjoint(analysis_type='aerothermoelastic') # 'aeroelastic') # 'aerothermoelastic') # 'aerothermal')
print('Created Object')

#dp.total_derivative()

optProb = Optimization("Stiffened Panel Aerothermoelastic Optimization", dp.objFunc)

optProb.addVarGroup("xvars", 112, "c", lower=0.001*np.ones(112), upper=0.1*np.ones(112), value=0.01)
optProb.addConGroup("con", 1, lower=33.4, upper=33.4)
optProb.addObj("obj")

print(optProb)

optOptions = {"IPRINT": -1}
opt = SLSQP(options=optOptions)
sol = opt(optProb, sens=dp.objGrad, storeHistory='opt_history.hst')

comm = MPI.COMM_WORLD
if comm.rank == 0:
    print(sol)
