#!/usr/bin/env python3
"""ATS input converter from version 1.4 to master branch"""

import sys, os
import numpy as np
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


def rename_element(xml, old, new):
    try:
        child = asearch.child_by_name(xml, old)
    except aerrors.MissingXMLError:
        pass
    else:
        child.setName(new)


def mafic_to_cap_pres(xml):
    """Converts mafic potential --> capillary pressure"""
    lc_list = asearch.find_path(xml, ["state", "initial conditions", "land cover types"], no_skip=True)
    for lc in lc_list:
        rename_element(lc, "mafic potential at fully closed stomata [Pa]",
                       "capillary pressure at fully closed stomata [Pa]")
        rename_element(lc, "mafic potential at fully open stomata [Pa]",
                       "capillary pressure at fully open stomata [Pa]")
            

def rooting_depth_fraction(xml):
    """rooting depth fraction --> root fraction"""
    try:
        rf = asearch.find_path(xml, ["state", "evaluators", "rooting_depth_fraction"], no_skip=True)
    except aerrors.MissingXMLError:
        pass
    else:
        etype = asearch.child_by_name(rf, "evaluator type")
        if etype.getValue() == "rooting depth fraction":
            etype.setValue("root fraction")
        rf.setName("root_fraction")


def has_surface_temp(sT):
    if asearch.child_by_name(sT, "evaluator type").getValue() == "independent variable":
        try:
            func_comp = asearch.find_path(sT, ["function", "function", "function-composition", "function1", "gradient"])
            print('found grad:')
            print(func_comp)
        except aerrors.MissingXMLError:
            return False
        else:
            val = func_comp.getValue()
            print('value = ', val)
            if np.allclose(np.array(val), np.array([1.0,0.0,0.0]), 1.e-6):
                # fix it!
                return True
            return False


def surface_temp(xml):
    """Many input files have a list where the surface temp was incorrectly set"""
    try:
        sT = asearch.find_path(xml, ["state", "evaluators", "surface-temperature"], no_skip=True)
    except aerrors.MissingXMLError:
        return

    if has_surface_temp(sT):
        for func in asearch.findall_path(sT, ["function", "function", "function-composition"]):
            f1 = asearch.child_by_name(func, "function1")
            f2 = asearch.child_by_name(func, "function2")
            if list(f1)[0].getName() == "function-linear" and \
               list(f2)[0].getName() != "function-linear":
                f2.setName("function1")
                f1.setName("function2")
                x0 = asearch.child_by_name(list(f1)[0], "x0")
                x0.setValue([86400., 0., 0.])
                grad = asearch.child_by_name(list(f1)[0], "gradient")
                grad.setValue([1.,0.,0.])
                y0 = asearch.child_by_name(list(f1)[0], "y0")
                y0.setValue(0.)

    try:
        snowmelt = asearch.find_path(xml, ["state", "evaluators", "snow-melt"])
    except aerrors.MissingXMLError:
        pass
    else:
        smr = asearch.child_by_name(snowmelt, "air-snow temperature difference [C]")

        try:
            snow_exp_temp = asearch.find_path(xml, ["state", "evaluators", "snow-expected_temperature"])
        except aerrors.MissingXMLError:
            # make one
            eval_list = asearch.find_path(xml, ["state", "evaluators"])
            snow_exp_temp = eval_list.sublist("snow-expected_temperature")
            snow_exp_temp.setParameter("evaluator type", 'string', "additive evaluator")
            snow_exp_temp.setParameter("dependencies", 'Array(string)', ["surface-air_temperature"])
            snow_exp_temp.setParameter("constant shift", 'double', -smr.getValue())
        else:
            pass

        try:
            snow_temp = asearch.find_path(xml, ["state", "evaluators", "snow-temperature"])
        except aerrors.MissingXMLError:
            # make one
            eval_list = asearch.find_path(xml, ["state", "evaluators"])
            snow_temp = eval_list.sublist("snow-temperature")
            snow_temp.setParameter("evaluator type", 'string', "secondary variable from function")
            snow_temp.setParameter("dependencies", 'Array(string)', ["snow-expected_temperature"])
            st_func = snow_temp.sublist("function")
            st_func.setParameter("function type", 'string', "composition")
            f1 = st_func.sublist("function1")
            f1.setParameter("function type", 'string', "linear")
            f1.setParameter("gradient", 'Array(double)', [1.0])
            f1.setParameter("y0", 'double', 273.15)

            f2 = st_func.sublist("function2")
            f2.setParameter("function type", 'string', "standard math")
            f2.setParameter("operator", 'string', "negative")
            f2.setParameter("shift", 'double', 273.15)
        else:
            pass
        snowmelt.remove(smr)


    try:
        can_temp = asearch.find_path(xml, ["state", "evaluators", "canopy-temperature"])
    except aerrors.MissingXMLError:
        # make one
        eval_list = asearch.find_path(xml, ["state", "evaluators"])
        can_temp = eval_list.sublist("canopy-temperature")
        can_temp.setParameter("evaluator type", 'string', "additive evaluator")
        can_temp.setParameter("dependencies", 'Array(string)', ["surface-air_temperature"])
        can_temp.setParameter("constant shift", 'double', -1.5)

    # lastly, clean up references to these new variables
    for p in asearch.findall_path(xml, ["snow temperature key",]):
        if p.getValue() == "surface-temperature":
            p.setValue("snow-temperature")
    for p in asearch.findall_path(xml, ["canopy temperature key",]):
        if p.getValue() == "surface-temperature":
            p.setValue("canopy-temperature")

    # also snow-evap needs new snow tmep
    try:
        snow_evap = asearch.find_path(xml, ["state", "evaluators", "snow-evaporation"])
    except aerrors.MissingXMLError:
        pass
    else:
        try:
            st = asearch.child_by_name(snow_evap, "surface temperature key")
        except aerrors.MissingXMLError:
            snow_evap.setParameter("surface temperature key", "snow-temperature")
        else:
            if st.getValue() == "surface-temperature":
                st.setValue("snow-temperature")
    

def transpiration_relperm(xml):
    try:
        trans = asearch.find_path(xml, ["state", "evaluators", "transpiration"])
    except aerrors.MissingXMLError:
        return

    tt = asearch.child_by_name(trans, "evaluator type")
    if tt.getValue() == "transpiration distribution via rooting depth":
        tt.setValue("transpiration distribution via relative permeability")

    try:
        wp = asearch.find_path(xml, ["state", "evaluators", "plant_wilting_factor"])
    except aerrors.MissingXMLError:
        pass
    else:
        asearch.find_path(xml, ["state", "evaluators"]).pop(wp)

            
def update(xml, transpiration_distribution=False):
    """generic update calls all needed things"""
    mafic_to_cap_pres(xml)
    rooting_depth_fraction(xml)
    surface_temp(xml)

    if transpiration_distribution:
        transpiration_relperm(xml)
    

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("infile", help="input filename")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-i", "--inplace", action="store_true", help="fix file in place")
    group.add_argument("-o", "--outfile", help="output filename")

    parser.add_argument("-t", "--transpiration-distribution", action="store_true",
                        help="Update to the relperm based transpiration distribution")
    
    args = parser.parse_args()

    # check for orig file
    print("Converting file: %s"%args.infile)
    xml = aio.fromFile(args.infile, True)
    update(xml, args.transpiration_distribution)
    if args.inplace:
        aio.toFile(xml, args.infile)
    else:
        aio.toFile(xml, args.outfile)
    sys.exit(0)
    

    
