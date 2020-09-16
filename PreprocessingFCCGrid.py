import os, sys, shutil
from dolfin import *

def RefineNearBoundary(mesh, nrefine):
    for j in range(nrefine):
      print("Refining ...")
      mesh.init(mesh.topology().dim()-1, mesh.topology().dim()) # Initialise facet to cell connectivity
      markers = MeshFunction("bool", mesh, mesh.topology().dim())
      for c in cells(mesh):
        for f in facets(c):
          if f.exterior():
            markers[c] = True
            break
      mesh = refine(mesh, markers)   
    print("Refined meshes: %d cells, %d vertices"%(mesh.num_cells(), mesh.num_vertices()))
    return mesh


comm = MPI.comm_world
nprocs = comm.Get_size()
if (nprocs>1):
    print('Support serial computation only!')
    sys.exit()


exists = os.path.isfile('DmriFemLib.py')
isupdate = False
if (exists==False or isupdate==True):
    if isupdate==True:
        os.system("rm DmriFemLib.py")
    print("Load pre-defined functions from GitHub")
    os.system("wget --quiet https://raw.githubusercontent.com/van-dang/DMRI-FEM-Cloud/master/DmriFemLib.py")
from DmriFemLib import *


"""# Working on the mesh"""
mesh_file= 'fcc_grid_v3'
file_dir= 'https://raw.githubusercontent.com/van-dang/DMRI-FEM-Cloud/mesh/'+mesh_file+'.msh.zip'
os.system('wget --quiet https://raw.githubusercontent.com/rochishnu00/DMRI/'+mesh_file+'.msh.zip')
zip_exists = os.path.isfile(mesh_file+".msh.zip")
mesh_file_exists = os.path.isfile(mesh_file)
if (zip_exists==False):
    os.system("wget "+file_dir)
if (mesh_file_exists==False):
    os.system("unzip -q "+mesh_file+".msh.zip")
os.system("dolfin-convert "+mesh_file+".msh "+mesh_file+".xml")

mesh_name = "fcc_grid_v3"
mymesh = Mesh(mesh_name+".xml");
GetPartitionMarkers(mesh_name+".msh", "pmk_"+mesh_name+".xml")
partition_marker = MeshFunction("size_t", mymesh, mymesh.topology().dim())
File("pmk_"+mesh_name+".xml")>>partition_marker
phase, partion_list = CreatePhaseFunc(mymesh, [], [], partition_marker)

CheckAndCorrectPeriodicity(mymesh, 0, 1e-6)
CheckAndCorrectPeriodicity(mymesh, 1, 1e-6)

print("Save phase function")
File("phase.pvd")<<phase

print("Partition markers:", partion_list)

print('Mesh file: ', mesh_file)

mesh = Mesh(mesh_file+".xml")

V_DG = FunctionSpace(mesh, 'DG', 0)

D0_array = [554.7e-6, 1664.2e-6, 554.7e-6]
IC_array = [1]
T2_array = [1e6]

# Variable tensor
dofmap_DG = V_DG.dofmap()
d00 = Function(V_DG); d01 = Function(V_DG); d02 = Function(V_DG)
d10 = Function(V_DG); d11 = Function(V_DG); d12 = Function(V_DG)
d20 = Function(V_DG); d21 = Function(V_DG); d22 = Function(V_DG)
T2 = Function(V_DG); disc_ic = Function(V_DG);
        

print('Setting parameters to %d cells'%(mesh.num_cells()))

T2.vector()[:]      = T2_array[0];
disc_ic.vector()[:] = IC_array[0];
d00.vector()[:]     = D0_array[0];
d11.vector()[:]     = D0_array[1];
d22.vector()[:]     = D0_array[2];

kappa_tensor = [[0,0,0],[0,0,0],[0,0,0]] 
kappa_tensor[0][1] = kappa_tensor[1][0] = 0.1e-3
kappa_tensor[1][2] = kappa_tensor[2][1] = 0.1e-3
CR = FunctionSpace(mymesh, "CR",1)
dofmap= CR.dofmap()
kappa = Function(CR)
D= mymesh.topology().dim()
mymesh.init(D-1,D) # Build connectivity between facets and cells
facet_f = MeshFunction("double", mymesh, D - 1)
for facet in facets(mymesh): 
    adjacent_cells = facet.entities(D) 
    if len(adjacent_cells)>1:
        domain_marker1 = partition_marker[adjacent_cells[0]]                            
        domain_marker2 = partition_marker[adjacent_cells[1]] 
        if domain_marker1!=domain_marker2:
            facet_f[facet] = kappa_tensor[domain_marker1][domain_marker2]  
    file = File('faces.pvd')
    file << facet_f     
    for cell in cells(mymesh):                                                         
        dofs = dofmap.cell_dofs(cell.index())                                          
        for i, facet in enumerate(facets(cell)):                                       
            kappa.vector()[dofs[i]] = facet_f[facet] 
    file = File('kappa.pvd')
    file << kappa
                             
'''
for cell in cells(mesh):
      cell_dof = dofmap_DG.cell_dofs(cell.index())
      T2.vector()[cell_dof]      = T2_array[0];
      disc_ic.vector()[cell_dof] = IC_array[0];
      d00.vector()[cell_dof]     = D0_array[0];
      d11.vector()[cell_dof]     = D0_array[0];
      d22.vector()[cell_dof]     = D0_array[0];
'''

ofile = 'files.h5';

for i in range(0, len(sys.argv)):
      arg = sys.argv[i];
      if arg=='-o':
            ofile = sys.argv[i+1];

filename, file_extension = os.path.splitext(ofile)


ofile = filename+'.h5'

print("Write to ",ofile)
f = HDF5File(mesh.mpi_comm(), ofile, 'w')
f.write(mesh, 'mesh')
f.write(T2, 'T2');  f.write(disc_ic, 'ic'); 
f.write(d00, 'd00'); f.write(d01, 'd01'); f.write(d02, 'd02')
f.write(d10, 'd10'); f.write(d11, 'd11'); f.write(d12, 'd12')
f.write(d20, 'd20'); f.write(d21, 'd21'); f.write(d22, 'd22')
f.write(kappa, 'kappa')
print("Done")
