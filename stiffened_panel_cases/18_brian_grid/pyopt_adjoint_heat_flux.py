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
        self.num_tacs_dvs = 112

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

    # There are 8 plate segments and 4 stiffener segments.
    # Since the overall panel is symmetric, they can be set equal to the opposite of one another
    # to reduce design variables
    def symmetryIndex(self, xInput):
        totalIndex = np.array([[0, 1, 2, 3, 3, 2, 1, 0, 4, 5, 5, 4], 
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]])

        xOutput = np.zeros(len(totalIndex[0]))

        for i in range(0, len(xInput)):
            for j in range(0, len(xOutput)):
                if i == totalIndex[0, j]:
                    xOutput[totalIndex[1, j]] = xInput[i]

        return xOutput


    def revSymmetryIndex(self, xInput):
        totalIndex = np.array([[0, 1, 2, 3, 3, 2, 1, 0, 4, 5, 5, 4], 
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]])

        xOutput = np.zeros(6)

        for i in range(0, len(xInput)):
            for j in range(0, len(totalIndex[0])):
                if i == totalIndex[1, j]:
                    xOutput[totalIndex[0, j]] += xInput[i]

        return xOutput

    # CAPS group assigments are all jumbled up. This maps them properly
    # First 8 indexes represent plate segments perpendicular to stiffeners
    # Last 4 indexes represent stiffener segments
    def designIndex(self, xInput):
        plateIndex = np.array([[0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7], 
        [16,23,40,47,64,71,88,95,103,111,15,22,39,46,63,70,87,94,102,110,14,21,38,45,62,69,86,93,101,109,13,20,37,44,61,68,85,92,100,108,12,19,36,43,60,67,84,91,99,107,11,18,35,42,59,66,83,90,98,106,10,17,34,41,58,65,82,89,97,105,8,9,32,33,56,57,80,81,96,104]])

        stiffenerIndex = np.array([[8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11],
        [72,73,74,75,76,77,78,79,48,49,50,51,52,53,54,55,24,25,26,27,28,29,30,31,0,1,2,3,4,5,6,7]])

        totalIndex = np.hstack((plateIndex, stiffenerIndex))
        xOutput = np.zeros(len(totalIndex[0]))

        for i in range(0, len(xInput)):
            for j in range(0, len(xOutput)):
                if i == totalIndex[0, j]:
                    xOutput[totalIndex[1, j]] = xInput[i]

        return xOutput


    def revDesignIndex(self, xInput):
        plateIndex = np.array([[0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7], 
        [16,23,40,47,64,71,88,95,103,111,15,22,39,46,63,70,87,94,102,110,14,21,38,45,62,69,86,93,101,109,13,20,37,44,61,68,85,92,100,108,12,19,36,43,60,67,84,91,99,107,11,18,35,42,59,66,83,90,98,106,10,17,34,41,58,65,82,89,97,105,8,9,32,33,56,57,80,81,96,104]])

        stiffenerIndex = np.array([[8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11],
        [72,73,74,75,76,77,78,79,48,49,50,51,52,53,54,55,24,25,26,27,28,29,30,31,0,1,2,3,4,5,6,7]])

        totalIndex = np.hstack((plateIndex, stiffenerIndex))
        xOutput = np.zeros( totalIndex[0, len(totalIndex[0]) - 1] + 1)

        for i in range(0, len(xInput)):
            for j in range(0, len(totalIndex[0])):
                if i == totalIndex[1, j]:
                    xOutput[totalIndex[0, j]] += xInput[i]
                    
        return xOutput

    def revdesignVarIndex(self, xInput):
        designVarIndex = np.array([111,16,110,102,94,87,70,63,46,39,22,103,15,109,101,93,86,69,62,45,38,21,95,14,108,100,92,85,68,61,44,37,20,88,13,107,99,91,84,67,60,43,36,19,71,12,106,98,90,83,66,59,42,35,18,64,11,105,97,89,82,65,58,41,34,17,47,10,104,96,80,81,56,57,32,33,8,40,9,23,79,53,29,5,76,52,28,4,75,51,27,55,3,74,50,26,2,73,49,25,1,72,31,48,24,0,7,78,54,30,6,77])
        xOutput = np.zeros(len(xInput))
        for i in range(0, len(xInput)):
            xOutput[designVarIndex[i]] = xInput[i]
        return xOutput

    def _build_model(self):

        thickness = 0.0001

        # Build the model
        model = FUNtoFEMmodel('wedge')
        plate = Body('plate', analysis_type=self.analysis_type, group=0,boundary=1)

        for i in range(self.num_tacs_dvs):
            plate.add_variable('structural',Variable('thickness '+ str(i),value=thickness,lower = 0.0001, upper = 0.01))

        #plate.add_variable('structural',Variable('thickness',value=thickness,lower = 0.01, upper = 0.1))
        model.add_body(plate)

        steady = Scenario('steady', group=0, steps=30)
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

    def objFunc(self, xdict):

        tInput1 = xdict["xvars"]
        tInput2 = self.symmetryIndex(tInput1)
        tInput3= self.designIndex(tInput2)

        self.model.set_variables(tInput3)
        self.driver.solve_forward()
        functions = self.model.get_functions()

        funcs = {}
        funcs["obj"] = functions[0].value
        funcs["con"] = functions[1].value

        # if self.comm.rank == 0:
        print('\n---------- FUNCTION SOLVE ----------')
        print('\nDesign Vars:       ', tInput1) 
        print('\nObjective Value:   ', functions[0].value)
        print('\nConstraint Value:  ', functions[1].value)

        fail = False

        return funcs, fail

    def objGrad(self, xdict, funcs): 

        tInput1 = xdict["xvars"]
        tInput2 = self.symmetryIndex(tInput1)
        tInput3= self.designIndex(tInput2)

        self.model.set_variables(tInput3)
        self.driver.solve_adjoint()
        grads = self.model.get_function_gradients()

        grad1_1 = np.array(grads[0][:])
        grad2_1 = np.array(grads[1][:])

        grad1_2 = self.revdesignVarIndex(grad1_1)
        grad1_3 = self.revDesignIndex(grad1_2)
        grad1_4 = self.revSymmetryIndex(grad1_3)

        grad2_2 = self.revdesignVarIndex(grad2_1)
        grad2_3 = self.revDesignIndex(grad2_2)
        grad2_4 = self.revSymmetryIndex(grad2_3)

        # if self.comm.rank == 0:
        print('\n---------- GRADIENT SOLVE ----------')
        print('\nDesign Vars:            ', tInput1) 
        print('\nObjective Grad Value:   ', grad1_4)
        print('\nConstraint Grad Value:  ', grad2_4)

        sens = {}
        sens = {
            "obj": {
                "xvars": [grad1_4]
            },
            "con": {
                "xvars": [grad2_4]
            },
        }

        fail = False

        return sens, fail


#==================================================================================================#

dp = wedge_adjoint(analysis_type='aerothermoelastic') # 'aeroelastic') # 'aerothermoelastic') # 'aerothermal')
print('Created Object')



optProb = Optimization("Stiffened Panel Aerothermoelastic Optimization", dp.objFunc)

optProb.addVarGroup("xvars", 6, "c", lower=0.00001*np.ones(6), upper=0.001*np.ones(6), value=0.0001)
# optProb.addVarGroup("xvars", 6, "c", lower=0.01*np.ones(6), upper=1*np.ones(6), value=0.1)

optProb.addConGroup("con", 1, lower=3, upper=4)
# optProb.addConGroup("con", 1, lower=300, upper=400)
optProb.addObj("obj")

comm = MPI.COMM_WORLD
# if comm.rank == 0:
print(optProb)

optOptions = {"IPRINT": -1}
opt = SLSQP(options=optOptions)
sol = opt(optProb, sens=dp.objGrad)

# if comm.rank == 0:
print(sol)
print('sol.xStar:  ', sol.xStar)