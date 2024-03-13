#!/usr/bin/env python3
"""ATS input converter from main development branch to tpetra development branch"""

import sys, os
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

            
def ics_to_constants(xml):
    try:
        ics = asearch.find_path(xml, ["state", "initial conditions"])
    except aerrors.MissingXMLError:
        pass
    else:
        ics.setName("constants")

def viscosity(xml):
    try:
        visc = asearch.find_path(xml, ["state", "evaluators", "viscosity_liquid"], no_skip=True)
    except aerrors.MissingXMLError:
        pass
    else:
        visc.setName("viscosity")


def inverse(xml):
    for inverse in asearch.findall_path(xml, ["PKs", "inverse"]):
        if inverse.isElement("preconditioning method"):
            pc = inverse.getElement("preconditioning method")
            if pc.getValue() == "boomer amg":
                pc.setValue("hypre: boomer amg")
                pcp = inverse.sublist("boomer amg parameters")
                pcp.setName("hypre: boomer amg parameters")


def change_eval_type(xml, name, old_type, new_type):
    try:
        e = asearch.find_path(xml, ["state", "evaluators", name], no_skip=True)
    except aerrors.MissingXMLError:
        return False
    else:
        etype = e.getElement("evaluator type")
        if etype.getValue() == old_type:
            etype.setValue(new_type)
            return True
        return False

def change_all_eval_type(xml, old_type, new_type):
    for eval_list in asearch.find_path(xml, ["state", "evaluators"], no_skip=True):
        if eval_list.getParameter("evaluator type").getValue() == old_type:
            eval_list.getParameter("evaluator type").setValue(new_type)
    
def wrm(xml):
    change_eval_type(xml, "saturation_liquid", "WRM", "wrm van Genuchten")
    change_eval_type(xml, "saturation_liquid", "water retention model", "wrm van Genuchten")
    change_eval_type(xml, "saturation_liquid", "wrm van Genuchten by material", "wrm van Genuchten")
    change_eval_type(xml, "saturation_gas", "WRM", "wrm van Genuchten")
    change_eval_type(xml, "saturation_gas", "water retention model", "wrm van Genuchten")
    change_eval_type(xml, "saturation_gas", "wrm van Genuchten by material", "wrm van Genuchten")
    change_eval_type(xml, "relative_permeability", "relative permeability, van Genuchten",
                     "relative permeability van Genuchten")
    change_eval_type(xml, "relative_permeability", "relative permeability, water retention model",
                     "relative permeability van Genuchten")
    change_eval_type(xml, "relative_permeability", "relative permeability van Genuchten by material",
                     "relative permeability van Genuchten")
    did_poro = change_eval_type(xml, "porosity", "compressible porosity", "compressible porosity linear by material")
    if did_poro:
        poro_eval = asearch.find_path(xml, ["state", "evaluators", "porosity"])
        if poro_eval.isElement("compressible porosity model parameters"):
            poro_eval.getElement("compressible porosity model parameters").setName("model parameters")

def hydraulic_conductivity(xml):
    try:
        rp = asearch.find_path(xml, ["state", "evaluators", "relative_permeability"], no_skip=True)
    except aerrors.MissingXMLError:
        return

    try:
        hc = asearch.find_path(xml, ["state", "evaluators", "relative_hydraulic_conductivity"], no_skip=True)
    except aerrors.MissingXMLError:
        hc = asearch.find_path(xml, ["state", "evaluators"], no_skip=True).sublist("relative_hydraulic_conductivity")
    else:
        return

    hc.setParameter("evaluator type", "string", "relative hydraulic conductivity")
    if rp.hasParameter("use surface rel perm"):
        hc.append(rp.pop("use surface rel perm"))
    if rp.hasParameter("density key"):
        hc.append(rp.pop("density key"))
    if rp.hasParameter("viscosity key"):
        hc.append(rp.pop("viscosity key"))

            
def perm_to_tensor(xml):
    try:
        e = asearch.find_path(xml, ["state", "evaluators", "permeability"], no_skip=True)
    except aerrors.MissingXMLError:
        pass
    else:
        etype = e.getParameter("evaluator type")
        if etype.getValue() == "independent variable":
            etype.setValue("independent variable tensor")
            e.setParameter("tensor rank", "int", 1)
        elif etype.getValue() == "independent variable constant":
            val = e.pop("value").getValue()
            etype.setValue("independent variable tensor")
            e.setParameter("constant in time", "bool", True)
            e.setParameter("tensor rank", "int", 1)
            dlist = e.sublist("function").sublist("domain")
            dlist.setParameter("region", "string", "computational domain")
            dlist.setParameter("component", "string", "cell")
            flist = dlist.sublist("function").sublist("function-constant")
            flist.setParameter("value", "double", val)



def bcs_with_units(xml):
    for bp in asearch.findall_path(xml, ["PKs", "boundary conditions", "boundary pressure"]):
        bp.setName("boundary pressure [Pa]")
    for bf in asearch.findall_path(xml, ["PKs", "boundary conditions", "outward water flux"]):
        bf.setName("outward water flux [mol m^-2 s^-1]")


def move_lc(xml):
    # move "snow melt rate" parameter to model parameters sublist
    eval_list = asearch.find_path(xml, ["state", "evaluators"], no_skip=True)
    for e in eval_list:
        if e.getParameter("evaluator type").getValue() == "snow melt rate":
            if e.hasParameter("snow melt rate [mm day^-1 C^-1]"):
                smr = e.pop("snow melt rate [mm day^-1 C^-1]")
                e.sublist("model parameters").append(smr)

    # move incident shortwave radiation parameters to model parameters
    for e in eval_list:
        if e.getParameter("evaluator type").getValue() == "incident shortwave radiation":
            if e.isElement("incident shortwave radiation parameters") and \
               not e.isElement("model parameters"):
                e.sublist("incident shortwave radiation parameters").setName("model parameters")
    
                
    # move land cover pars to model pars list
    consts = asearch.find_path(xml, ["state", "constants"])
    try: 
        lc = consts.pop("land cover types")
    except aerrors.MissingXMLError:
        pass
    else:
        lc.setName("land cover parameters")
        mp = asearch.find_path(xml, ["state", "model parameters"], no_skip=True)
        mp.append(lc)

    # add land cover pars to model pars list in evalutors
    for e in eval_list:
        if e.getParameter("evaluator type").getValue() in ["potential evapotranspiration, Priestley-Taylor",
                                                           "radiation balance, surface, snow, and canopy",
                                                           "area fractions, two components",
                                                           "area fractions, three components",
                                                           "albedos, two components",
                                                           "albedos, three components",
                                                           "transpiration distribution, relative permeability",
                                                           "root fraction",]:
            if not e.hasParameter("model parameters"):
                e.setParameter("model parameters", "string", "land cover parameters")

        elif e.getParameter("evaluator type").getValue() == "snow melt rate":
            if not e.hasParameter("additional model parameters"):
                e.setParameter("additional model parameters", "Array(string)", ["land cover parameters",])


def by_material(xml):
    for e in asearch.find_path(xml, ["state", "evaluators"], no_skip=True):
        if e.getParameter("evaluator type").getValue().endswith("by material"):
            etype = e.getParameter("evaluator type")
            etype.setValue(etype.getValue().replace("by material", "by region"))
                

def update(xml):
    """generic update calls all needed things""" 
    ics_to_constants(xml)
    viscosity(xml)
    inverse(xml)
    wrm(xml)
    perm_to_tensor(xml)
    bcs_with_units(xml)
    hydraulic_conductivity(xml)

    change_all_eval_type(xml, "additive evaluator", "additive")
    change_all_eval_type(xml, "multiplicative evaluator", "multiplicative")
    change_all_eval_type(xml, "radiation balance, surface and canopy", "radiation balance, surface, snow, and canopy")
    change_all_eval_type(xml, "subgrid albedos, two components", "albedos, two components")
    change_all_eval_type(xml, "subgrid albedos, three components", "albedos, three components")

    move_lc(xml)
    by_material(xml)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("infile", help="input filename")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-i", "--inplace", action="store_true", help="fix file in place")
    group.add_argument("-o", "--outfile", help="output filename")

    args = parser.parse_args()

    print("Converting file: %s"%args.infile)
    xml = aio.fromFile(args.infile, True)
    update(xml)
    if args.inplace:
        aio.toFile(xml, args.infile)
    else:
        aio.toFile(xml, args.outfile)
    sys.exit(0)
    

    
