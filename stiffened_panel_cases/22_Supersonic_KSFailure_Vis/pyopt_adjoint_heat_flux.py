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

    def __init__(self, analysis_type='aerothermoelastic'):
        print('start')

        self.analysis_type = analysis_type

        # cruise conditions
        self.v_inf = 171.5                  # freestream velocity [m/s]
        self.rho = 0.01841                  # freestream density [kg/m^3]
        self.cruise_q = 12092.5527126               # dynamic pressure [N/m^2]
        self.grav = 9.81                            # gravity acc. [m/s^2]
        self.thermal_scale = 0.5 * self.rho * (self.v_inf)**3

        self.maximum_mass = 40.0 
        self.num_tacs_dvs = 12

        # Set up the communicators
        n_tacs_procs = 1

        comm = MPI.COMM_WORLD
        self.comm = comm

        world_rank = comm.Get_rank()
        if world_rank < n_tacs_procs:
            color = 55
            key = world_rank
        else:
            color = MPI.UNDEFINED
            key = world_rank
        self.tacs_comm = comm.Split(color,key)
        # Set up the FUNtoFEM model for the TOGW problem
        self._build_model()
        self.ndv = len(self.model.get_variables())
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

        self.obj_scale = 0.0106
        self.con_scale = 3.35

        self.optHist = open('optHist.txt', 'w')
        self.optHistAll = open('optHistAll.txt', 'w')
        


    def _build_model(self):

        thickness = 0.001

        # Build the model
        model = FUNtoFEMmodel('wedge')
        plate = Body('plate', analysis_type=self.analysis_type, group=0,boundary=1)

        for i in range(self.num_tacs_dvs):
            plate.add_variable('structural',Variable('thickness '+ str(i),value=thickness,lower = 0.0001, upper = 0.01))

        model.add_body(plate)

        steady = Scenario('steady', group=0, steps=75)
        temp = Function('ksfailure',analysis_type='structural') # Objective
        steady.add_function(temp)

        mass = Function('mass',analysis_type='structural',adjoint=False) # Constraint
        steady.add_function(mass)

        model.add_scenario(steady)

        self.model = model


    def objFunc(self, xdict):

        tInput1 = xdict["xvars"]

        self.model.set_variables(tInput1)
        self.driver.solve_forward()
        functions = self.model.get_functions()

        func1 = functions[0].value * 1.0/self.obj_scale
        func2 = functions[1].value * 1.0/self.con_scale

        funcs = {}
        funcs["obj"] = func1
        funcs["con"] = func2

        fail = False

        if self.comm.rank == 0:
            print('\n---------- FUNCTION SOLVE ----------', flush=True)
            print('\nDesign Vars:       ', tInput1, flush=True) 
            print('\nObjective Value:   ', func1, flush=True)
            print('\nConstraint Value:  ', func2, flush=True)

            print('\n---------- FUNCTION SOLVE ----------', flush=True, file=self.optHist)
            print('\nDesign Vars:       ', tInput1, flush=True, file=self.optHist) 
            print('\nObjective Value:   ', func1, flush=True, file=self.optHist)
            print('\nConstraint Value:  ', func2, flush=True, file=self.optHist)

        print('\nself.comm.rank:      ', self.comm.rank, flush=True, file=self.optHistAll)
        print('\n---------- FUNCTION SOLVE ----------', flush=True, file=self.optHistAll)
        print('\nDesign Vars:       ', tInput1, flush=True, file=self.optHistAll) 
        print('\nObjective Value:   ', func1, flush=True, file=self.optHistAll)
        print('\nConstraint Value:  ', func2, flush=True, file=self.optHistAll)

        return funcs, fail


    def objGrad(self, xdict, funcs): 

        tInput1 = xdict["xvars"]

        self.model.set_variables(tInput1)
        self.driver.solve_adjoint()
        grads = self.model.get_function_gradients()

        grad1 = np.array(grads[0][:]) * 1.0/self.obj_scale
        grad2 = np.array(grads[1][:]) * 1.0/self.con_scale

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

        if self.comm.rank == 0:
            print('\n---------- GRADIENT SOLVE ----------', flush=True)
            print('\nDesign Vars:            ', tInput1, flush=True) 
            print('\nObjective Grad Value:   ', grad1, flush=True)
            print('\nConstraint Grad Value:  ', grad2, flush=True)

            print('\n---------- GRADIENT SOLVE ----------', flush=True, file=self.optHist)
            print('\nDesign Vars:            ', tInput1, flush=True, file=self.optHist) 
            print('\nObjective Grad Value:   ', grad1, flush=True, file=self.optHist)
            print('\nConstraint Grad Value:  ', grad2, flush=True, file=self.optHist)
        
        print('\nself.comm.rank:           ', self.comm.rank, flush=True, file=self.optHistAll)
        print('\n---------- GRADIENT SOLVE ----------', flush=True, file=self.optHistAll)
        print('\nDesign Vars:            ', tInput1, flush=True, file=self.optHistAll) 
        print('\nObjective Grad Value:   ', grad1, flush=True, file=self.optHistAll)
        print('\nConstraint Grad Value:  ', grad2, flush=True, file=self.optHistAll)

        return sens, fail


#==================================================================================================#

dp = wedge_adjoint(analysis_type='aerothermoelastic') # 'aeroelastic') # 'aerothermoelastic') # 'aerothermal')
print('Created Object')



optProb = Optimization("Stiffened Panel Aerothermoelastic Optimization", dp.objFunc)

optProb.addVarGroup("xvars", 12, "c", lower=0.0001*np.ones(12), upper=0.01*np.ones(12), value=0.001)

optProb.addConGroup("con", 1, lower=1, upper=1)
optProb.addObj("obj")

comm = MPI.COMM_WORLD
# if comm.rank == 0:
print(optProb)

optOptions = {"IPRINT": -1}
opt = SLSQP(options=optOptions)
sol = opt(optProb, sens=dp.objGrad)

if comm.rank == 0:
    print(sol)
    print(sol, file=dp.optHist)
    print('\nsol.xStar:  ', sol.xStar)
    print('\nsol.xStar:  ', sol.xStar, file=dp.optHist)
print(sol, file=dp.optHistAll)
print('\nsol.xStar:  ', sol.xStar, file=dp.optHistAll)

dp.optHist.close()
dp.optHistAll.close()