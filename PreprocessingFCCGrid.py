import os, sys, shutil
from dolfin import *

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
mesh_name= "fcc_grid_v3"

os.system('wget --quiet https://raw.githubusercontent.com/rochishnu00/DMRI/'+mesh_name+'.geo')

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


mesh_name = "fcc_grid_v3"
    
# Create mesh from geo file by gmsh
os.system('gmsh -3 '+mesh_name+'.geo -o '+mesh_name+'.msh')

# Convert .msh to .xml using dolfin-convert
os.system('dolfin-convert '+mesh_name+'.msh '+mesh_name+'.xml')

#os.system('wget --quiet https://raw.githubusercontent.com/rochishnu00/DMRI/'+mesh_name+'.xml')

mymesh = Mesh(mesh_name+".xml");  

GetPartitionMarkers(mesh_name+".msh", "pmk_"+mesh_name+".xml")

partition_marker = MeshFunction("size_t", mymesh, mymesh.topology().dim())

File("pmk_"+mesh_name+".xml")>>partition_marker
    
phase, partion_list = CreatePhaseFunc(mymesh, [], [], partition_marker)
    
File("Phase.pvd")<<phase

print("Save phase function")
File("phase.pvd")<<phase

print("Partition markers:", partion_list)

ofile = 'files.h5';

for i in range(0, len(sys.argv)):
      arg = sys.argv[i];
      if arg=='-o':
            ofile = sys.argv[i+1];

filename, file_extension = os.path.splitext(ofile)

ofile = filename+'.h5'

f = HDF5File(mymesh.mpi_comm(), ofile, 'w')
f.write(mymesh, 'mesh');  f.write(T2, 'T2'); f.write(disc_ic, 'ic'); f.write(phase, 'phase');
print("Write to ", ofile)
