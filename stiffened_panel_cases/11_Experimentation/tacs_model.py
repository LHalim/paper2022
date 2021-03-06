from __future__ import print_function
import os

from pprint import pprint
import numpy as np
from mpi4py import MPI

from tacs import TACS, functions, constitutive, elements, pyTACS, problems
from pyfuntofem.tacs_interface import TacsSteadyInterface

class wedgeTACS(TacsSteadyInterface):
    def __init__(self, comm, tacs_comm, model, n_tacs_procs):
        super(wedgeTACS,self).__init__(comm, tacs_comm, model)

        assembler = None
        self.tacs_proc = False
        if comm.Get_rank() < n_tacs_procs:
            self.tacs_proc = True

            # Instantiate FEASolver
            structOptions = {
                'printtiming':True,
            }

            bdfFile = os.path.join(os.path.dirname(__file__), 'stiffPanel4Thermal.dat')
            FEASolver = pyTACS(bdfFile, options=structOptions, comm=tacs_comm)

            # Material properties
            rho = 2780.0        # density kg/m^3
            E = 73.1e9          # Young's modulus (Pa)
            nu = 0.33           # Poisson's ratio
            kcorr = 5.0/6.0     # shear correction factor
            ys = 324.0e6        # yield stress
            specific_heat = 920.096
            cte = 24.0e-6
            kappa = 230.0

            # tInput = 0.001*np.ones(112)

            softPanelT = 0.001*np.ones(6)
            tInput2 = self.symmetryIndex(softPanelT)
            tInput3= self.designIndex(tInput2)

            def elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):
                elemIndex = kwargs['propID'] - 1
                # t = tInput[elemIndex]
                t = tInput3[elemIndex]

                prop = constitutive.MaterialProperties(rho=rho, specific_heat=specific_heat,
                                                       E=E, nu=nu, ys=ys, cte=cte, kappa=kappa)
                con = constitutive.IsoShellConstitutive(prop, t=t, tNum=dvNum)

                refAxis = np.array([1.0, 0.0, 0.0])

                elemList = []
                transform = None
                for elemDescript in elemDescripts:
                    if elemDescript in ['CQUAD4', 'CQUADR']:
                        elem = elements.Quad4ThermalShell(transform, con)
                    else:
                        print("Uh oh, '%s' not recognized" % (elemDescript))
                    elemList.append(elem)

                # Add scale for thickness dv
                scale = [100.0]
                return elemList, scale

            # # Set up elements and TACS assembler
            # FEASolver.initialize(elemCallBack)
            # assembler = FEASolver.assembler

            # # Add back pressure
            # forces = assembler.createVec()
            # force_array = forces.getArray()

            # # Panel pressure loading
            # force_array[2::7] += 0.005*np.cos(np.radians(5))
            # force_array[0::7] += -0.005*np.sin(np.radians(5))
            
            # assembler.setBCs(forces)

        self._initialize_variables(assembler, thermal_index=6)

        self.initialize(model.scenarios[0],model.bodies)
    
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



    def post_export_f5(self):
        flag = (TACS.OUTPUT_CONNECTIVITY |
                TACS.OUTPUT_NODES |
                TACS.OUTPUT_DISPLACEMENTS |
                TACS.OUTPUT_STRAINS |
                TACS.OUTPUT_STRESSES |
                TACS.OUTPUT_EXTRAS)
        f5 = TACS.ToFH5(self.assembler, TACS.BEAM_OR_SHELL_ELEMENT, flag)
        f5.writeToFile('stiffPanel.f5')
