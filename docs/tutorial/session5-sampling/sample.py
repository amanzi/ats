# On EES-16 servers, make sure that an ats module is loaded to access these modules.
# Otherwise, add the tools/utils folder in your ats source directory to your PYTHONPATH environment variable
import atsxml
import ats_xdmf
# On EES-16 servers, load the matk module: module load matk
# Otherwise, directions to obtain and install matk are at http://matk.lanl.gov/installation.html
from matk import matk

# The model that MATK will provide parameters to as a dictionay in the argument (pars)
def model(pars,hostname,processor):

    # Modify base ats xml input file with pars dictionary and run ats
    m = atsxml.get_root('../test7-v_fwd.xml')
    atsxml.replace_by_path(m,['base_porosity','rest domain','value'],pars['poro_m'])
    atsxml.replace_by_path(m,['base_porosity','peat','value'],pars['poro_p'])
    atsxml.replace_by_path(m,['permeability','rest domain','value'],10**pars['perm_m'])
    atsxml.replace_by_path(m,['permeability','peat','value'],10**pars['perm_p'])
    atsxml.run(m,nproc=4,stdout='stdout.out',stderr='stdout.err',cpuset=processor) 

    # Read results from ats visualization files
    vf = ats_xdmf.VisFile()
    # Collect point at middle of polygon 1 m deep
    # x,z = 7.17946807, 4.65764252
    # index for this location is 1733, see below how to find this
    # Create output dictionary that matches MATK observations
    out = {}
    out['Sl'] = vf.get('saturation_liquid', vf.cycles[-1])[1733]
    out['T'] = vf.get('temperature', vf.cycles[-1])[1733]
    vf.close()

    # Return simulated values of interest
    return out

# Create host dictionary so that cpu sets can be explicitly defined for ats runs
# The host (dictionary key) 'dum' will be ignored in this case
# The list of strings will be used as the cpu sets for the runs
hosts = {'dum':['0,1,2,3','4,5,6,7','8,9,10,11']}
# Create MATK object specifying the 'model' function above as the model
p = matk(model=model)

# Add some parameters
# Mineral soil porosity
p.add_par('poro_m',min=0.586, max=0.606, value=0.596)
# Peat porosity
p.add_par('poro_p', min=0.866, max=0.886, value=0.876)
# Mineral soil permeability
p.add_par('perm_m',min=-13.5, max=-12.5, value=-13)
# Peat permeability
p.add_par('perm_p', min=-12.5, max=-11.5, value=-12)

# Create observations
# Saturation
p.add_obs('Sl',value=0.5)
# Temperature
p.add_obs('T',value=273.)

# Create parameter study where nvals indicats the number of values for each parameter
# Parameter values will be spaced evenly over parameter ranges
s = p.parstudy(nvals=[2,2,2,3])
# Save parameter combinations in file
s.savetxt('sample.txt')
# Run sampleset in individual folders with base name 'run'
# logfile will contain results of sampling in progress, including error messages for failed runs
# outfile will be written at the end with samples in correct order
s.run(cpus=hosts,workdir_base='run',outfile='sample.out',logfile='sample.log')

## How to determine index for cell in ATS mesh
#import mesh
#import numpy
## cd to directory with completed ATS run
## Collect cell centroid coordinates, reads visdump_mesh.h5 in current directory by default
#ecs = mesh.meshElemCentroids()
## Look at unique values of x along transect
#numpy.unique(ecs[:,0])
## Collect coordinates in center of polygon (x = 7.17946807)
#x7 = numpy.where(abs(ecs[:,0]-7.17946807)<0.0001)
## Find max z at center of polygon
#numpy.max(ecs[x7][:,2])
## Find zs at 1 meter deep
#z1 =numpy.where(abs(ecs[:,2]-4.6597)<0.01)
## Find index at intersection of xs at center of polygon and zs around 4.6597 m deep
#obs = numpy.intersect1d(x7,z1)[0]
## Use value of obs (1733) to collect results in "model" function above




