import sys, os
import subprocess
from ConfigParser import SafeConfigParser as config_parser
import netCDF4
import unittest

    
_filenames = {"extrude_one_layer":"Mesh3D_OneLayer",
              "extrude_uniform":"Mesh3D_2mSoil",
              "extrude_homogeneous_uniform":"Mesh3D_Homogeneous2mSoil",
              "extrude_variable":"Mesh3D_VariableSoil",
              "extrude_homogeneous_variable":"Mesh3D_HomogeneousVariableSoil"}


def runExe(dirname, exes=None):
    cwd = os.getcwd()
    if exes is None:
        exes = _filenames.keys()

    try:
        os.chdir(dirname)
        for exe in exes:
            exo = _filenames[exe]+".exo"
            try:
                os.remove(exo)
            except:
                pass
        
            log = _filenames[exe]+".log"

            executable = os.path.join("..", "..", exe)
            print "Running: %s in %s"%(executable, os.getcwd())
            
            with open(log, 'w') as stdout:
                subprocess.call(executable, stdout=stdout)

    finally:
        os.chdir(cwd)



def loadConfig(dirname):
    cp = config_parser()
    cp.read(os.path.join(dirname, "test.cfg"))
    return cp

def openMesh(meshname):
    return netCDF4.Dataset(meshname, 'r')


def findDirectories():
    dirs = os.listdir('.')
    dirs = [d for d in dirs if 'run_test.py' not in d]
    return dirs



# --------------------------------------------------------------------------------
# actual tests

class Mesh(unittest.TestCase):
    def setUp(self):
        self.fid = openMesh(os.path.join(self.dirname, ".".join([self.meshname,"exo"])))
    def tearDown(self):
        self.fid.close()

def generateTestDimension(attr, value):
    def test(self):
        self.assertEqual(len(self.fid.dimensions[attr]), int(value))
    return test

def generateMeshClass(dirname, meshname, config):
    class MyMesh(Mesh):
        pass

    setattr(MyMesh, "dirname", dirname)
    setattr(MyMesh, "meshname", meshname)

    sec_name = ".".join([meshname, "dimensions"])
    for dname in config.options(sec_name):
        setattr(MyMesh, "test_%s"%dname,
                generateTestDimension(dname, config.get(sec_name, dname)))
    MyMesh.__name__ = dirname+"_"+meshname
    MyMesh.__module__ = ""
    return MyMesh

def generateSuites(dirname):
    suites = []
    config = loadConfig(dirname)
    sections = config.sections()
    meshes = set(sec.split(".")[0] for sec in sections)

    for meshname in meshes:
        suite = unittest.TestLoader().loadTestsFromTestCase(generateMeshClass(dirname, meshname, config))
        suites.append(suite)
    return suites
    


if __name__ == "__main__":
    suites = []
    for d in findDirectories():
        runExe(d)
        suites.append(unittest.TestSuite(generateSuites(d)))
    suite = unittest.TestSuite(suites)
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
        
    
