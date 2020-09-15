# -*- coding: utf-8 -*-
# This program solves the Bloch-Torrey equation applied to computational diffusion MRI using 
# the finite element method coupled with the theta-method for the spatial discretization.

# The scope of usage: 
# (1) Single domains, Multilayered structures, manifolds
# (2) Membrane permeability for internal interfaces
#     Artificial permeability at the external interfaces
# (3) pure homogeneous Neumann BCs, (4) pseudo-periodic BCs

# Copyright (C) 2019 Van-Dang Nguyen (vdnguyen@kth.se)

# This file is part of DOLFIN.

# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.

# First added:  2017-10-10
# Last changed: 2019-04-25

# This demo is maintained by Van-Dang Nguyen
# Please report possible problems to vdnguyen@kth.se

# Disable warnings
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger("rothemain.rothe_utils")
logging.getLogger('UFL').setLevel(logging.WARNING)
logging.getLogger('FFC').setLevel(logging.WARNING)


import warnings
warnings.filterwarnings("ignore")

import os, sys, shutil

from dolfin import *
parameters['ghost_mode'] = 'shared_facet' 
parameters['allow_extrapolation'] = True



## default parameters
g0, g1, g2 = 0, 1, 0; kcoeff = 3e-3; porder = 1; Nsteps = 100; bvalue = 1000; kappa = 1e-5;
delta, Delta = 10600, 43100; T2 = 1e16; k = 200; nskip = 5; is_input_dt = 0; is_input_b = 0; is_input_q = 0;
is_kcoeff_from_file = 1; is_T2_from_file = 1; is_IC_from_file = 0; IsDomainPeriodic = False; IsDomainMultiple = False
PeriodicDir = [0, 0, 0]; phase = None; ic = None
## end default parameters

comm = MPI.comm_world
rank = comm.Get_rank()

try:
    ## Input parameters from command lines
        for i in range(0, len(sys.argv)):
                arg = sys.argv[i];
                if arg=='-f':
                        ffile = sys.argv[i+1];
                        if rank==0:
                            print('input file:', ffile)
                if arg=='-M':
                        IsDomainMultiple = int(sys.argv[i+1]);
                        if rank==0:
                            print('IsDomainMultiple:', IsDomainMultiple)
                if arg=='-IsPeriodic' or arg=='-isperiodic':
                        IsDomainPeriodic = int(sys.argv[i+1]);
                        if rank==0:
                            print('IsDomainPeriodic:', IsDomainPeriodic)
                if arg=='-N':
                        Nsteps = int(sys.argv[i+1]);
                        if rank==0:
                            print('Nsteps:', Nsteps)
                if arg=='-b':
                        is_input_b = 1;
                        bvalue = float(sys.argv[i+1]);
                        if rank==0:
                            print('bvalue:', bvalue)
                if arg=='-q':
                        is_input_q = 1;
                        qvalue = float(sys.argv[i+1]);
                        if rank==0:
                            print('qvalue:', bvalue)
                if arg=='-p':
                        kappa = float(sys.argv[i+1]);
                        if rank==0:
                            print('permeability:', kappa)
                if arg=='-D':                       
                        Delta = float(sys.argv[i+1]);
                        if rank==0:
                            print('Delta:', Delta)
                if arg=='-d':
                        delta = float(sys.argv[i+1]);
                        if rank==0:
                            print('delta:', delta)
                if arg=='-K':
                        is_kcoeff_from_file = 0
                        if rank==0:
                            print("Reading diffusion coefficient from command line")
                        kcoeff = float(sys.argv[i+1]);
                        if rank==0:
                            print('diffusion coefficient:', kcoeff)
                if arg=='-k':
                        is_input_dt = 1;
                        k = float(sys.argv[i+1]);
                        if rank==0:
                            print('time step size:', k)
                if arg=='-T2':
                        is_T2_from_file = 0;
                        T2 = float(sys.argv[i+1]);
                        if rank==0:
                            print('T2: ', T2)
                if arg=='-gdir':
                        g0 = float(sys.argv[i+1]);
                        g1 = float(sys.argv[i+2]);
                        g2 = float(sys.argv[i+3]);
                        if rank==0:
                            print('(g0, g1, g2):',g0, g1, g2)
                if arg=='-pdir':
                        PeriodicDir[0] = int(sys.argv[i+1]);
                        PeriodicDir[1] = int(sys.argv[i+2]);
                        PeriodicDir[2] = int(sys.argv[i+3]);
                        if rank==0:
                            print('PeriodicDir=[', PeriodicDir[0], PeriodicDir[1], PeriodicDir[2],']')                        
except:
        print('Something goes wrong with the inputs!')

"""#Load pre-defined functions"""

exists = os.path.isfile('DmriFemLib.py')
isupdate = False
if (exists==False or isupdate==True):
    if isupdate==True:
        os.system("rm DmriFemLib.py")
    if rank==0:    
        print("Load pre-defined functions from GitHub")
    os.system("wget --quiet https://raw.githubusercontent.com/van-dang/DMRI-FEM-Cloud/master/DmriFemLib.py")

from DmriFemLib import *
"""# Solve the Bloch-Torrey equation"""

mymesh = Mesh()
myf = HDF5File(mymesh.mpi_comm(),ffile, 'r')
myf.read(mymesh, 'mesh', False)

CheckAndCorrectPeriodicity(mymesh, 0, 1e-6)
CheckAndCorrectPeriodicity(mymesh, 1, 1e-6)

V_DG = FunctionSpace(mymesh, 'DG', 0)
dofmap_DG = V_DG.dofmap()

mri_simu = MRI_simulation()
mri_para = MRI_parameters()
#bvalues = [500]
bvalues = [0]
for bvalue in bvalues:
    #################################################################################
    #########################  Pre-defined parameters ###############################
    mri_para.stype='PGSE'                                # sequence type
    mri_para.bvalue = bvalue;                            # bvalue
    mri_para.delta, mri_para.Delta = 10000, 10000           # time sequence
    mri_para.set_gradient_dir(mymesh, 1, 1, 0)           # gradient direction

    mri_para.T = mri_para.Delta+mri_para.delta
    mri_para.fs_sym = sp.Piecewise(
                    (  1., mri_para.s < mri_para.delta ),
                    (  0., mri_para.s < mri_para.Delta ),
                    ( -1., mri_para.s < mri_para.T ),
                    (  0., True )  
                ) 

    mri_para.Apply()
    mri_simu.k = 100;                                    # time-step size
    mri_simu.nskip = 5;                                  # frequency to print ouputs
    mydomain = MyDomain(mymesh, mri_para)
    mydomain.phase = phase 
    mydomain.PeriodicDir = [1, 1, 0];             # Direction of the periodicity
    mydomain.IsDomainPeriodic = True            # Confirm if the mesh if periodic
    mydomain.IsDomainMultiple = True              # Confirm if the mesh is multiple
    ################################################################################
    # mydomain.kappa_e_scalar = 5e-4;
    # mydomain.kappa = 1;
    kappa_tensor = [[0,0,0],[0,0,0],[0,0,0]] 
    kappa_tensor[0][1] = kappa_tensor[1][0] = 0.1e-3
    kappa_tensor[1][2] = kappa_tensor[2][1] = 0.1e-3
    CR = FunctionSpace(mydomain.mymesh, "CR",1)
    dofmap= CR.dofmap()
    kappa = Function(CR)
    D= mydomain.mymesh.topology().dim()
    mydomain.mymesh.init(D-1,D) # Build connectivity between facets and cells
    facet_f = MeshFunction("double", mydomain.mymesh, D - 1)
    for facet in facets(mydomain.mymesh): 
      adjacent_cells = facet.entities(D) 
      if len(adjacent_cells)>1:
        domain_marker1 = partition_marker[adjacent_cells[0]]                            
        domain_marker2 = partition_marker[adjacent_cells[1]] 
        if domain_marker1!=domain_marker2:
          facet_f[facet] = kappa_tensor[domain_marker1][domain_marker2]  
    file = File('faces.pvd')
    file << facet_f     
    for cell in cells(mydomain.mymesh):                                                         
      dofs = dofmap.cell_dofs(cell.index())                                          
      for i, facet in enumerate(facets(cell)):                                       
        kappa.vector()[dofs[i]] = facet_f[facet] 
    file = File('kappa.pvd')
    file << kappa
    mydomain.kappa=avg(kappa)
    mydomain.Apply()
    # Impose the diffusion coefficient
    # mydomain.D  = 3e-3;
    DO_array=[554.7e-6, 1664.2e-6, 554.7e-6]
    #Variable diffusion tensor
    V_DG=mydomain.V_DG; dofmap_DG=V_DG.dofmap();
    d00= Function(V_DG); d01=Function(V_DG); d02=Function(V_DG);
    d10= Function(V_DG); d11=Function(V_DG); d12=Function(V_DG);
    d20= Function(V_DG); d21=Function(V_DG); d22=Function(V_DG);
    mydomain.ImposeDiffusionTensor(d00,d01,d02,d10,d11,d12,d20,d21,d22);

    #################################################################################

    # linsolver = PETScLUSolver("mumps")
    linsolver = KrylovSolver("bicgstab")
    linsolver.parameters["absolute_tolerance"] = 1e-4
    linsolver.parameters["relative_tolerance"] = 1e-4
    linsolver.parameters["maximum_iterations"] = 10000
    
    mri_simu.solve(mydomain, mri_para, linsolver)

    ctext = ''; # text to comment in the log file.
    PostProcessing(mydomain, mri_para, mri_simu, plt, ctext)
