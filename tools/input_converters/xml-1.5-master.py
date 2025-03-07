#!/usr/bin/env python3
"""ATS input converter from version 1.5 to master branch"""

import sys, os
import numpy as np
import copy
try:
    amanzi_xml = os.path.join(os.environ["AMANZI_SRC_DIR"], "tools","amanzi_xml")
except KeyError:
    pass
else:
    if amanzi_xml not in sys.path:
        sys.path.append(amanzi_xml)

from amanzi_xml.utils import search as asearch
from amanzi_xml.utils import io as aio
from amanzi_xml.utils import errors as aerrors
from amanzi_xml.common import parameter, parameter_list


def enforceDtHistory(xml):
    """Find and revert the timestep from file option, moving it to the cycle driver list."""
    ti_file_pars = None
    for ti in asearch.findall_name(xml, "time integrator"):
        if asearch.child_by_name(ti, "timestep controller type").getValue() == "from file":
            for child in ti:
                if child.getName().startswith("timestep controller") and \
                   child.getName().endswith("parameters") and \
                   child.getName() != "timestep controller from file parameters":
                    new_type = child.getName()[len("timestep controller "):-len(" parameters")]
                    asearch.child_by_name(ti, "timestep controller type").setValue(new_type)
                    break
            ti_file_pars = ti.pop("timestep controller from file parameters")

    if ti_file_pars is not None:
        # push the filename into the cycle driver list instead.  We
        # just use the last one found -- they should all be the same!
        cycle_driver_list = asearch.child_by_name(xml, "cycle driver")
        tsm_list = cycle_driver_list.sublist("timestep manager")
        tsm_list.setParameter("prescribed timesteps file name", "string", asearch.child_by_name(ti_file_pars, "file name").getValue())

def timeStep(xml):
    """Many parameters changed "time step" --> "timestep" """
    def fix(xml, pname):
        for p in asearch.findall_name(xml, pname):
            p.setName(pname.replace("time step", "timestep"))

    fix(xml, "time step reduction factor")
    fix(xml, "time step control method")
    fix(xml, "time step increase factor")        
    fix(xml, "time step cut factor")        
    fix(xml, "time step cut threshold")        
    fix(xml, "time step increase threshold")        
    fix(xml, "max time step")        
    fix(xml, "min time step")        
    fix(xml, "max time step [s]")        
    fix(xml, "min time step [s]")        
    fix(xml, "max time step (s)")        
    fix(xml, "min time step (s)")        
    fix(xml, "initial time step")
    fix(xml, "initial time step [s]")        
    fix(xml, "initial time step (s)")        
    fix(xml, "max valid change in saturation in a time step [-]")
    fix(xml, "max valid change in ice saturation in a time step [-]")
    fix(xml, "subcycling target time step [s]")
    fix(xml, "time step")

    for ti in asearch.findall_name(xml, "timestep controller fixed parameters"):
        if ti.isElement("initial timestep [s]"):
            ti.getElement("initial timestep [s]").setName("timestep [s]")

    for ti_type in ["standard", "smarter"]:
        for ti in asearch.findall_name(xml, f"timestep controller {ti_type} parameters"):
            if not ti.isElement("initial timestep [s]"):
                ti.setParameter("initial timestep [s]", "double", 1.0)


def tensorPerm(xml):
    for eval_list in asearch.find_path(xml, ["state", "evaluators"], True):
        ename = eval_list.getName()
        if ename == "permeability" or ename.endswith("-permeability"):
            tensorPerm_(eval_list)

def tensorPerm_(perm):
    # set the type to tensor
    ptype = perm.getElement("evaluator type")
    if ptype.getValue() == "independent variable":
        ptype.setValue("independent variable tensor")
    elif ptype.getValue() == "independent variable constant":
        ptype.setValue("independent variable tensor")

        # create the function list
        value = perm.getElement("value").getValue()
        flist = perm.sublist("function").sublist("domain")
        flist.setParameter("region", "string", "computational domain")
        flist.setParameter("component", "string", "cell")
        flist.sublist("function").sublist("function-constant").setParameter("value", 'double', value)

    else:
        return
            
    # set the rank
    if not perm.isElement("tensor rank"):
        if perm.isElement("permeability type"):
            tt = perm.getElement("permeability type")
            tt.setName("tensor type")
            if tt.getValue() == "full tensor":
                tt.setValue("full symmetric")
            elif tt.getValue() == "diagonal tensor":
                tt.setValue("diagonal")
        else:
            perm.setParameter("tensor type", "string", "scalar")



def initialConditionsList(xml):
    for pk in xml.sublist("PKs"):
        if pk.isElement("initial condition"):
            pk.sublist("initial condition").setName("initial conditions")
            
            
def update(xml):
    #enforceDtHistory(xml)
    timeStep(xml)
    tensorPerm(xml)
    initialConditionsList(xml)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("infile", help="input filename")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-i", "--inplace", action="store_true", help="fix file in place")
    group.add_argument("-o", "--outfile", help="output filename")
    args = parser.parse_args()

    # check for orig file
    print("Converting file: %s"%args.infile)
    xml = aio.fromFile(args.infile, True)
    update(xml)
    if args.inplace:
        aio.toFile(xml, args.infile)
    else:
        aio.toFile(xml, args.outfile)
    sys.exit(0)
