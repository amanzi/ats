import sys,os
sys.path.append(os.path.join(os.environ['ATS_DIR'],'tools','meshing_ats'))
import meshing_ats

import numpy as np


# set up the surface mesh
xmax = 100.0
ymax = 1.0
slope = 0.05
nsurf_cells = 100

x_c = np.linspace(0.0, xmax, nsurf_cells + 1)
y_c = np.array([0.0, ymax])




topsoil_max_thickness = 0.1
topsoil_nmesh_layers = 2
restsoil_zones_thickness = [0.4, 0.6, 2.0]
restsoil_nmesh_layers = [4, 3, 2]


Xc, Yc = np.meshgrid(x_c, y_c)
Xc = Xc.flatten()
Yc = Yc.flatten()

restsoil_thickness = np.sum(restsoil_zones_thickness)
zmax_restsoil = restsoil_thickness + slope*xmax

# connectivity
conn = []
for i in range(nsurf_cells):
    conn.append([i, i+1, nsurf_cells + i + 2, nsurf_cells + i + 1])

#top soil thickness function
def topsoil_thickness_fun(x,y):
    return topsoil_max_thickness - \
           np.minimum(2*topsoil_max_thickness*np.sin(np.pi*x/xmax)**2, topsoil_max_thickness)

Zc = topsoil_max_thickness + restsoil_thickness - slope*Xc - \
     np.minimum(2*topsoil_max_thickness*np.sin(np.pi*Xc/xmax)**2, topsoil_max_thickness)

coords = np.array([Xc, Yc, Zc])

#make 2D mesh
m2 = meshing_ats.Mesh2D(coords.transpose(), conn)

types = ['function',] + ['constant',]*len(restsoil_zones_thickness)
thicknesses = [topsoil_thickness_fun,] + restsoil_zones_thickness
ncells = [topsoil_nmesh_layers,] + restsoil_nmesh_layers
mat_ids = [10000,] + [20000,]*(len(restsoil_zones_thickness) - 1) + [30000,]

#make 3D mesh
m3 = meshing_ats.Mesh3D.extruded_Mesh2D(m2, types, thicknesses, ncells, mat_ids)
#write 3D mesh
m3.write_exodus("meshing_ats_example.exo")
