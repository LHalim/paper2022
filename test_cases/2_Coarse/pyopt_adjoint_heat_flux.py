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
#from pyfuntofem.massoud_body import *

from tacs_model import wedgeTACS
#from pyOpt import Optimization,SLSQP
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

        self.var_scale        = np.ones(self.ndv,dtype=TransferScheme.dtype)
        self.struct_tacs = solvers['structural'].assembler

    def _build_model(self):

        thickness = 0.015
        # Build the model
        model = FUNtoFEMmodel('wedge')
        plate = Body('plate', analysis_type=self.analysis_type, group=0,boundary=1)
        plate.add_variable('structural',Variable('thickness',value=thickness,lower = 0.01, upper = 0.1))
        model.add_body(plate)

        steady = Scenario('steady', group=0, steps=100)
        steady.set_variable('aerodynamic',name='AOA',value=0.0,lower=-15.0,upper=15.0)
        temp = Function('ksfailure',analysis_type='structural') #temperature
        steady.add_function(temp)

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

        return 

    def verification_test(self,epsilon=1e-6,steps=100):
        
        steady = self.model.scenarios[0]
        bodies = self.model.bodies
        body = self.model.bodies[0]
                
        fail = self.driver.solve_forward()
        fail = self.driver.solve_adjoint()

        #solve_forward() 1
        self.driver._distribute_variables(steady, bodies)
        self.driver._distribute_functions(steady, bodies)
        self.driver._initialize_forward(steady, bodies)
        self.driver._update_transfer()
        #self.driver._solve_steady_forward(steady)
        for step in range(1,steps+1):
            fail = self.driver.solvers['flow'].iterate(steady, bodies, step)
        self.driver._post_forward(steady, bodies)
       
        # Store Output
        if body.transfer is not None:
            # Aeroelastic Terms
            body.aero_loads_copy = body.aero_loads.copy()
        if body.thermal_transfer is not None:
            # Aerothermal Terms
            body.aero_heat_flux_mag_copy = body.aero_heat_flux_mag.copy() 

        #solve_adjoint()
        self.driver._distribute_variables(steady, bodies)
        self.driver._distribute_functions(steady, bodies)
        self.driver._initialize_adjoint_variables(steady, bodies)
        self.driver._initialize_adjoint(steady, bodies)

        if body.transfer is not None:
            # Aeroelastic Terms
            body.dLdfa = np.ones(shape=body.dLdfa.shape)       
        if body.thermal_transfer is not None:
            # Aerothermal Terms
            body.dQdfta = np.ones(shape=body.dQdfta.shape)

        #self.driver._solve_steady_adjoint(steady)
        for step in range(1,steps+1):
            fail = self.driver.solvers['flow'].iterate_adjoint(steady, bodies, step)
        self.driver._post_adjoint(steady, bodies)

        # Perturb and Get Adjoint Product
        if body.transfer is not None:
            # Aeroelastic Terms
            adjoint_product = 0.0
            #body.aero_disps_pert = np.random.uniform(size=body.aero_disps.shape)
            body.aero_disps_pert = np.ones(shape=body.aero_disps.shape)
            body.aero_disps_pert[::2] = 0.5
            body.aero_disps += epsilon*body.aero_disps_pert
            adjoint_product += np.dot(body.dGdua[:, 0], body.aero_disps_pert) 
        if body.thermal_transfer is not None:
            # Aerothermal Terms
            adjoint_product_t = 0.0
            #body.aero_temps_pert = np.random.uniform(size=body.aero_temps.shape)
            body.aero_temps_pert = np.ones(shape=body.aero_temps.shape)
            body.aero_temps_pert[::2] = 0.5
            body.aero_temps += epsilon*body.aero_temps_pert
            adjoint_product_t += np.dot(body.dAdta[:, 0], body.aero_temps_pert)
        
        #solve_forward() 2
        self.driver._distribute_variables(steady, bodies)
        self.driver._distribute_functions(steady, bodies)
        self.driver._initialize_forward(steady, bodies)
        self.driver._update_transfer()
        #self.driver._solve_steady_forward(steady)
        for step in range(1,steps+1):
            fail = self.driver.solvers['flow'].iterate(steady, bodies, step)
        self.driver._post_forward(steady, bodies)

        #Finite Difference
        if body.transfer is not None:
            # Aeroelastic Terms
            fd_product = 0.0
            fd = (body.aero_loads - body.aero_loads_copy)/epsilon
            fd_product += np.dot(fd, body.dLdfa[:, 0])
            print('FUN3D FUNtoFEM adjoint result (elastic):           ', adjoint_product)
            print('FUN3D FUNtoFEM finite-difference result (elastic): ', fd_product)
        if body.thermal_transfer is not None:
            #Thermal Terms
            fd_product_t = 0.0
            fd = (body.aero_heat_flux_mag - body.aero_heat_flux_mag_copy)/epsilon
            fd_product_t += np.dot(fd, body.dQdfta[:, 0])
            print('FUN3D FUNtoFEM adjoint result (thermal):           ', adjoint_product_t)
            print('FUN3D FUNtoFEM finite-difference result (thermal): ', fd_product_t)
            """
            print('body.aero_heat_flux = ', body.aero_heat_flux)
            print('body.aero_heat_flux_copy = ', body.aero_heat_flux_copy)
            print('fd = ', fd)
            print('fd_product_t = ', fd_product_t)
            print('body.dQdfta = ', body.dQdfta)
            """
        if body.transfer is not None and body.thermal_transfer is not None:    
            # Total
            print('FUN3D FUNtoFEM adjoint result (total):           ', (adjoint_product + adjoint_product_t))
            print('FUN3D FUNtoFEM finite-difference result (total): ', (fd_product + fd_product_t))        

        return 
        

################################################################################
dp = wedge_adjoint(analysis_type='aerothermal') # 'aeroelastic') # 'aerothermoelastic') # 'aerothermal')
print('created object')

# print('VERIFICATION TEST')
# Error = dp.verification_test(epsilon=1e-6)
# print('FINISHED VERIFICATION TEST')

dp.total_derivative()
