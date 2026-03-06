"""A script to strip out an ATS checkpoint file and write a "logically
structured" h5 version for use as ICs by PFLOTRAN
"""

import sys,os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'],'tools','utils'))

import numpy as np
import h5py
import mesh
import ats_xdmf

def ats_to_pflotran_ic_h5(filename, directory=".", output_filename="pflotran_ic.h5"):
    ixyz = mesh.meshCentroidsStructuredOrdering_3D(directory=directory)

    with h5py.File(os.path.join(directory, filename),'r') as fin:
        ic_pressure = fin['pressure.cell.0'][:][ixyz['id']]
        ic_temperature = fin['temperature.cell.0'][:][ixyz['id']]

        with h5py.File(os.path.join(directory, output_filename),'w') as fout:
            fout.create_dataset("pressure", data=ic_pressure)
            fout.create_dataset("temperature", data=ic_temperature)
            fout.create_dataset("x", data=ixyz['x'])
            fout.create_dataset("y", data=ixyz['y'])
            fout.create_dataset("z", data=ixyz['z'])

def ats_to_pflotran_bcs_h5(directory=".", output_filename="pflotran_bcs.h5"):
    ixy = mesh.meshCentroidsStructuredOrdering_3D(order=["x",], filename="visdump_surface_mesh.h5",
                                                  directory=directory)

    vf = ats_xdmf.VisFile(directory, filename="visdump_surface_data.h5",
                          output_time_unit='s')
    keys = vf.cycles
    times = vf.times
    dat = vf.d

    with h5py.File(os.path.join(directory, output_filename),'w') as fout:
        fout.create_dataset("time [s]", data=np.array(times))
        T = fout.create_group("surface temperature [K]")
        for i,k in enumerate(keys):
            T.create_dataset("%d"%i, data=dat['surface-temperature.cell.0'][k][:][ixy['id']])

        flx = fout.create_group("outward molar flux [mol m^-2 s^-1]")

        # need the face area
        face_areas = mesh.meshElemVolume(filename="visdump_surface_mesh.h5", directory=directory)

        for i,k in enumerate(keys):
            flux_dat = dat['surface_subsurface_flux.cell.0'][k][:]
            flux_dat = flux_dat / face_areas
            flx.create_dataset("%d"%i, data=flux_dat[ixy['id']])

    vf.close()


if __name__ == "__main__":
    checkp = sys.argv.pop(-1)
    if not (checkp.startswith("checkpoint") and checkp.endswith(".h5")):
        print "Usage: python ats-vis-to-structured-h5.py checkpointXXXXX.h5"
        sys.exit(-1)

    ats_to_pflotran_ic_h5(checkp)
    ats_to_pflotran_bcs_h5()
    sys.exit(0)

    
