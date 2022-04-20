import sys
sys.path.append("../")
import meshing_ats

m2 = meshing_ats.Mesh2D.read_VTK("lobster_full.vtk")

# layer thicknesses, except for the last entry which is the coordinate of the bottom
types = ['constant',]*5 + ['snapped',]
thicknesses = [0.2, 0.8, 1.0, 5.0, 13.0, -40.0]
ncells = [int(0.2/0.02), int(0.8 / 0.02), int(1./0.05), int(5./0.2), int(13/1.0), int(20/2.0)]
print "total ncells = ", sum(ncells)
mat_ids = [10000,]+[20000,]*5

m3 = meshing_ats.Mesh3D.extruded_Mesh2D(m2, types, thicknesses, ncells, mat_ids)
m3.write_exodus("lobster_full-one_block.exo", "one block")

# m3 = meshing_ats.Mesh3D.extruded_Mesh2D(m2, types, thicknesses, ncells, mat_ids)
# m3.write_exodus("lobster_full-not_duplicated.exo", "n blocks, not duplicated")

# m3 = meshing_ats.Mesh3D.extruded_Mesh2D(m2, types, thicknesses, ncells, mat_ids)
# m3.write_exodus("lobster_full-duplicated.exo", "n blocks, duplicated")

# m3 = meshing_ats.Mesh3D.extruded_Mesh2D(m2, types, thicknesses, ncells, mat_ids)
# m3.write_exodus("lobster_full-repeated.exo", "one block, repeated")


