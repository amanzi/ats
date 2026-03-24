"""Scripts much of the work to take a run and include it in the test suite.

Expects the following environmental variables to be defined:
  AMANZI_SRC_DIR = /path/to/amanzi
  ATS_SRC_DIR = /path/to/ats

There are a few subtleties to regression tests.  The biggest one is
the issue of adaptive timestepping.  Changing the timestep history
often results in a larger change in the solution (due to
non-convergence in time) than a reasonable epsilon to ensure success
of the regression test.  Even things like machine-precision roundoff
differences can grow sufficiently to result in a single step taking
one more or less nonlinear iterations, then driving the adaptive
timestepper to change the timestep size and therefore failing the
regression test.

This script does the following: it runs a new test, then parses the
logfile to generate a history of timestep sizes.  It then creates a
new input file which drives the same simulation with the same dT.
This simulation becomes the new test.
"""

from __future__ import print_function
from __future__ import division

import regression_tests
import argparse
import subprocess
import time
import textwrap
import sys
import os
import h5py

sys.path.append(os.path.join(os.environ['AMANZI_SRC_DIR'], 'tools', 'amanzi_xml'))
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], 'tools', 'utils'))
import parse_logfile
import ats_xdmf
from amanzi_xml.utils import errors as aerrors
from amanzi_xml.utils import io as aio
from amanzi_xml.utils import search as asearch


_txtwrap = textwrap.TextWrapper(width=78, subsequent_indent=4*" ")


def dirname(testname):
    return testname+".regression.orig"

def stdout(testname):
    return os.path.join(dirname(testname),testname+".stdout")


def commandline_options():
    """Process the command line arguments and return them as a dict."""
    parser = argparse.ArgumentParser(description='Set up and ATS regression test')

    parser.add_argument('-c', '--write_config', action='store_true',
                        help='write test to config file')
    parser.add_argument('-e', '--executable', nargs=1, default=None,
                        help='path to executable to use for testing')
    parser.add_argument('-m', '--mpiexec', nargs=1, default=None,
                        help="path to the executable for mpiexec (mpirun, etc)"
                        "on the current machine.")
    parser.add_argument('-n', '--num_processors', nargs=1, default=None, type=int,
                        help="number of processors on which to run this test")
    parser.add_argument('config_filename', metavar='CONFIG_FILENAME', type=str,
                        help="name of the configure file")
    parser.add_argument('testfiles', metavar='TEST_FILES', nargs='+', type=str,
                        help="list of tests (must all be in the same subdirectory)")

    options = parser.parse_args()
    return parser.parse_args()


def run(mpiexec, np, executable, testname, testlog):
    """Build up the run command and run the job as a subprocess."""
    command = []
    if np is not None and np is not 1:
        if not mpiexec:
            raise RuntimeError("Multiple process simulation requested but no valid mpiexec was provided.")
        command.append(mpiexec)
        command.append("-np")
        command.append(str(np))

    command.append(executable)
    command.append("--xml_file=../{0}.xml".format(testname))

    test_directory = os.getcwd()
    run_directory = os.path.join(test_directory, dirname(testname))
    os.mkdir(run_directory)
    os.chdir(run_directory)

    print("    cd {0}".format(run_directory), file=testlog)
    print("    {0}".format(" ".join(command)), file=testlog)
    with open(testname + ".stdout", 'w') as run_stdout:
        proc = subprocess.call(command,
                               stdout=run_stdout,
                               stderr=subprocess.STDOUT)

    if proc:
        raise RuntimeError("Unsuccesful run of the base run -- see {0}/{1}".format(test_directory,stdout(testname)))
    os.chdir(test_directory)
    
    
def build_dt_history(testname, testlog):
    """Reads the stdout of the run to get timestep history then writes them to hdf5 file"""
    with open(stdout(testname),'r') as fid:
        good, bad = parse_logfile.parse_file(fid)

    if not os.path.isdir('data'):
        os.mkdir('data')
        
    with h5py.File(os.path.join('data', "{0}_dts.h5".format(testname)), 'w') as fid:
        fid.create_dataset("timesteps", data=86400.*good[:,2]) # note conversion from days to seconds

    print("Wrote file: data/{0}_ats.h5".format(testname), file=testlog)


def find_vis_cycles(testname):
    """Parses the output of the run to collect the cycles that were written."""
    try:
        vf = ats_xdmf.VisFile(dirname(testname))
    except RuntimeError:
        try:
            vf = ats_xdmf.VisFile(dirname(testname), filename="visdump_surface_data.h5")
        except RuntimeError:
            raise RuntimeError("Unable to open {0}/visdump_data.h5 (or the surface equivalent)".format(testname))

    keys = vf.cycles
    vf.close()
    return [int(k) for k in keys]



def write_new_xml(testname, testlog):
    """Adds the needed magic to a new xml file."""
    print("Renaming: {0}.xml to {0}_orig.xml".format(testname), file=testlog)
    os.rename(testname+".xml", testname+"_orig.xml")

    xml = aio.fromFile(testname+"_orig.xml", True)

    # vis by list of cycles
    try:
        vis = asearch.getElementByNamePath(xml, "visualization")
    except aerrors.MissingXMLError:
        pass
    else:
        for i in range(len(vis)):
            vis.pop(vis[0].get('name'))
        vis.setParameter("cycles", "Array(int)", find_vis_cycles(testname))
        vis.setParameter("file name base", "string", "visdump")

    try: 
        viss = asearch.getElementByNamePath(xml, "visualization surface")
    except aerrors.MissingXMLError:
        pass
    else:
        for i in range(len(viss)):
            viss.pop(viss[0].get('name'))
        viss.setParameter("cycles", "Array(int)", find_vis_cycles(testname))
        viss.setParameter("file name base", "string", "visdump_surface")

    # update timestep controller, nonlinear solvers
    for ti in asearch.generateElementByNamePath(xml, "time integrator"):
        asearch.generateElementByNamePath(ti, "limit iterations").next().setValue(100)
        asearch.generateElementByNamePath(ti, "diverged tolerance").next().setValue(1.e10)

        asearch.getElementByNamePath(ti, "timestep controller type").setValue("from file")
        ts_hist = ti.sublist("timestep controller from file parameters")
        ts_hist.setParameter("file name", "string", "../data/{0}_dts.h5".format(testname))

    print("Writing: {0}.xml".format(testname), file=testlog)
    aio.toFile(xml, "{0}.xml".format(testname))



    
config_header = """
[suites]

[default-test-criteria]
default = 1.0e-12 absolute

"""
def write_to_config(testname, options, testlog):
    """Opens the config file, creating it if need be, and appends the test."""
    config_file = os.path.split(options.config_filename)[1]
    if not os.path.isfile(config_file):
        print("New config file : {0}".format(config_file), file=testlog)
        with open(config_file,'w') as fid:
            fid.write(config_header)

    with open(config_file, 'a') as fid:
        fid.write("[{0}]\n".format(testname))
        if (options.num_processors is not None and options.num_processors is not 1):
            fid.write("np = {0}\n".format(options.num_processors))
        fid.write("\n")
        
            

        
def generate_new_test(options, testlog):
    """The main function to generate new tests."""

    # set up executables
    mpiexec = regression_tests.check_for_mpiexec(options, testlog)
    executable = regression_tests.check_for_executable(options, testlog)

    for testfile in options.testfiles:
        # determine the test name
        testname = testfile[:]
        if testname.endswith('.xml'):
            testname = testname[:-4]

        # run and process    
        run(mpiexec, options.num_processors, executable, testname, testlog)
        build_dt_history(testname, testlog)
        write_new_xml(testname, testlog)
        if options.write_config:
            write_to_config(testname, options, testlog)
        
        
if __name__ == '__main__':
    options = commandline_options()
    print(options.config_filename)
    test_directory = os.path.dirname(options.config_filename)
    print("Generating tests in: {0}".format(test_directory))
    testlog = regression_tests.setup_testlog(_txtwrap)
    print("Generating tests in: {0}".format(test_directory), file=testlog)

    os.chdir(test_directory)
    generate_new_test(options, testlog)
    
    
    
        



