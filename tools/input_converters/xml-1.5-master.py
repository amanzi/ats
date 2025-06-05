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
import fix_chemistry_ts_control

def moveLC(xml):
    ics = asearch.find_path(xml, ["state", "initial conditions"])

    if ics.isElement("land cover types"):
        lc = ics.pop("land cover types")
        model_pars = xml.sublist("state").sublist("model parameters")
        model_pars.append(lc)

def viscosityRelationType(xml):
    for par in asearch.findall_path(xml, ["state", "evaluators", "viscosity relation type"]):
        par.setName("viscosity type")

def molarMass(xml):
    for par in asearch.findall_path(xml, ["state", "evaluators", "molar mass"]):
        par.setName("molar mass [kg mol^-1]")
        

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


def _pkTreeIter(tree, recurse=True):
    yield tree
    for el in tree:
        if el.getName() != "PK type":
            for pk in _pkTreeIter(el, recurse):
                yield pk
            

def coupledFlowTransportMPC(xml):
    pk_tree = asearch.find_path(xml, ["cycle driver", "PK tree"])
    pk_list = asearch.find_path(xml, ["PKs",])

    transport_pk_types = ["transport ATS", "surface subsurface transport"]
    flow_pk_types = ["richards flow", "richards steady state", "overland flow, pressure basis",
                     "subsurface permafrost", "icy surface", "permafrost model", "coupled water",
                     "operator split coupled water", "operator split permafrost",]
    
    for pk in _pkTreeIter(list(pk_tree)[0]):
        if pk.getElement("PK type").getValue() == "subcycling MPC":
            subpk_types = [el.getElement("PK type").getValue() for el in pk if el.getName() != "PK type"]
            print(f'Found a subcycling MPC with subpk types: {subpk_types}')

            if len(subpk_types) == 2 and \
               subpk_types[0] in flow_pk_types and \
               subpk_types[1] in transport_pk_types:
                print(f'... changing to type "coupled flow and transport"')
                pk.getElement("PK type").setValue("coupled flow and transport")
                pk_list.sublist(pk.getName()).getElement("PK type").setValue("coupled flow and transport")


def removeTc99(xml):
    # Tc-99 not a valid name
    for cnames in asearch.findall_path(xml, ["component names",]):
        old_names = cnames.getValue()
        new_names = []
        for name in cnames.getValue():
            if name == 'Tc-99':
                new_names.append('Tracer1')
            elif '-' in name:
                new_names.append(name.replace('-', ''))
            else:
                new_names.append(name)
        cnames.setValue(new_names)
    
                
def removeVerboseObject(xml):
    pm = asearch.parent_map(xml)
    
    # let global verbosity control for tests
    for vo in asearch.findall_path(xml, ["verbose object",]):
        pm[vo].remove(vo)


def fixTransportPK(pk, evals_list):
    if pk.isElement("domain name"):
        domain = pk.getElement("domain name").getValue()
    else:
        domain = "domain"
        
    # spatial discretization order --> advection spatial discretization order
    if not pk.isElement("advection spatial discretization order") and pk.isElement("spatial discretization order"):
        order = pk.getElement("spatial discretization order")
        order.setName("advection spatial discretization order")
    order = pk.getElement("advection spatial discretization order").getValue()
    if order == 1 and pk.isElement("reconstruction"):
        pk.pop("reconstruction")

    # moves and changes format of molecular diffusion coefficient
    if pk.isElement("molecular diffusion") and not pk.isElement("molecular diffusivity [m^2 s^-1]"):
        md = pk.pop("molecular diffusion")
        if md.isElement("aqueous names"):
            names = md.getElement("aqueous names").getValue()
            vals = md.getElement("aqueous values").getValue()
            print("found MD for names:", names)
            if len(names) > 0 and names[0] != '':
                diff = pk.sublist("molecular diffusivity [m^2 s^-1]")
                for name, val in zip(names, vals):
                    diff.setParameter(name.strip(), "double", val)

    # inverse list in diffusion list
    if pk.isElement("inverse") and pk.isElement("diffusion"):
        inv = pk.pop("inverse")
        pk.sublist("diffusion").append(inv)
    elif pk.isElement("inverse"):
        pk.pop("inverse")

    # internal subcycling is dead
    if pk.isElement("transport subcycling"):
        pk.pop("transport subcycling")

    # runtime diagnostics is dead
    if pk.isElement("runtime diagnostics: regions"):
        pk.pop("runtime diagnostics: regions")

    if pk.isElement("material properties"):
        mat_props = pk.pop("material properties").sublist("domain")

        # moves dispersion coefficient to MDM type
        if mat_props.isElement("model"):
            if domain == "domain":
                disp_key = "dispersion_coefficient"
            else:
                disp_key = domain + "-dispersion_coefficient"

            if not evals_list.isElement(disp_key):
                all_zero = True
                disp_list = evals_list.sublist(disp_key)
                disp_list.setParameter("evaluator type", "string", "dispersion tensor")
                disp_list2 = disp_list.sublist("mechanical dispersion parameters").sublist(domain)
                if domain == "domain":
                    disp_list2.setParameter("region", "string", "computational domain")
                else:
                    disp_list2.setParameter("region", "string", domain + " domain")

                disp_type = mat_props.getElement("model").getValue()
                if disp_type == "scalar":
                    disp_list2.setParameter("mechanical dispersion type", "string", "isotropic")
                    params_list = mat_props.sublist(f"parameters for {disp_type}")
                    params_list.setName("isotropic parameters")
                    disp_list2.append(params_list)

                    for par in params_list:
                        if par.getValue() != 0.0:
                            all_zero = False
                else:
                    disp_list2.setParameter("mechanical dispersion type", "string", disp_type)
                    params_list = mat_props.sublist(f"parameters for {disp_type}")
                    params_list.setName(f"{disp_type} parameters")
                    disp_list2.append(params_list)

                    for par in params_list:
                        if par.getValue() != 0.0:
                            all_zero = False

                if all_zero:
                    # all params are zero, remove dispersivity compleletely
                    evals_list.pop(disp_key)

        # moves/reformats tortuosity
        if mat_props.isElement("aqueous tortuosity"):
            pk.sublist("tortuosity [-]").setParameter("aqueous", "double", mat_props.getElement("aqueous tortuosity").getValue())
        if mat_props.isElement("gaseous tortuosity"):
            pk.sublist("tortuosity [-]").setParameter("gaseous", "double", mat_props.getElement("gaseous tortuosity").getValue())

    # component molar masses gets units, moves to sublist
    if pk.isElement("component molar masses"):
        mms = pk.pop("component molar masses")
        names = pk.getElement("component names").getValue()
        for name, mm in zip(names, mms):
            pk.sublist("molar mass [kg mol^-1]").setParameter(name, "double", mm)

    # look for a water_content, define if needed
    if domain == "domain":
        lwc_key = "water_content"
    else:
        lwc_key = f"{domain}-water_content"

    print(f'searching for LWC = {lwc_key}')
        
    if not evals_list.isElement(lwc_key):
        print(' ... is not eval')
        lwc_list = evals_list.sublist(lwc_key)
        lwc_list.setParameter("evaluator type", "string", "multiplicative evaluator")

        if "surface" in domain:
            print('adding surface quantities for LWC')
            lwc_list.setParameter("dependencies", "Array(string)", [f"{domain}-ponded_depth", f"{domain}-molar_density_liquid", f"{domain}-cell_volume"])
        else:            
            print('adding subsurface quantities for LWC')
            lwc_list.setParameter("dependencies", "Array(string)", ["saturation_liquid", "molar_density_liquid", "porosity", "cell_volume"])

    # boundary conditions->concentration --> bondary conditions->mole fraction
    bcs_list = pk.sublist('boundary conditions')
    if bcs_list.isElement('concentration') and not bcs_list.isElement('mole fraction'):
        bcs_list.sublist('concentration').setName('mole fraction')

    if bcs_list.isElement('mole fraction'):
        mf = bcs_list.sublist('mole fraction')
        for sublist in mf:
            if sublist.isElement('boundary concentration function') and not sublist.isElement('boundary mole fraction function'):
                sublist.sublist('boundary concentration function').setName('boundary mole fraction function')
    

            
    # remove dead keys
    def removeDead(keyname, default_names, remove_eval=False):
        if isinstance(default_names, str):
            default_names = [default_names,]
        keyname_key = keyname + " key"

        if pk.isElement(keyname_key):
            key = pk.getElement(keyname_key).getValue()
            if key in default_names:
                pk.pop(keyname_key)
                if remove_eval and evals_list.isElement(key):
                    evals_list.pop(key)

        keyname_key = keyname_key + " suffix"
        default_names = [dn.split('-')[-1] for dn in default_names]
        if pk.isElement(keyname_key):
            key = pk.getElement(keyname_key).getValue()
            if key in default_names:
                pk.pop(keyname_key)

    removeDead("porosity", ["porosity",])
    removeDead("porosity", ["surface-porosity", "surface-one"], True)
    removeDead("molar density liquid", ["surface-molar_density_liquid", "molar_density_liquid"])
    removeDead("saturation", ["surface-saturation_liquid", "surface-ponded_depth", "saturation_liquid"])
    removeDead("saturation liquid", ["surface-saturation_liquid", "surface-ponded_depth", "saturation_liquid"])
    removeDead("flux", ["surface-water_flux", "water_flux"])
    removeDead("water flux", ["surface-water_flux", "water_flux"])
    if pk.isElement("flux_key"):
        pk.pop("flux_key")
    if pk.isElement("molar_density_key"):
        pk.pop("molar_density_key")

    if pk.isElement("number of liquid components"):
        pk.pop("number of liquid components")
    if pk.isElement("number of gaseous components") and pk.getElement("number of gaseous components").getValue() == 0:
        pk.pop("number of gaseous components")
    if pk.isElement("number of aqueous components") and pk.isElement("component names") and \
       pk.getElement("number of aqueous components").getValue() == len(pk.getElement("component names").getValue()):
        pk.pop("number of aqueous components")
    if pk.isElement("PK origin"):
        pk.pop("PK origin")
    if pk.isElement("solver"):
        pk.pop("solver")

    if pk.isElement("physical models and assumptions"):
        pk.pop("physical models and assumptions")
    if pk.isElement("flow mode"):
        pk.pop("flow mode")


def fixTransportPKs(xml):
    evals_list = asearch.find_path(xml, ["state", "evaluators"], True)
    pk_tree = asearch.find_path(xml, ["cycle driver", "PK tree"])
    pk_list = asearch.find_path(xml, ["PKs",])

    for pk in _pkTreeIter(list(pk_tree)[0]):
        pk_type = pk.getElement("PK type")

        if pk_type.getValue() == "transport ATS":
            print("fixing transport pk")
            fixTransportPK(pk_list.sublist(pk.getName()), evals_list)
        elif pk_type.getValue() == "transport ats":
            pk_type.setValue("transport ATS")
            print("fixing transport pk")
            fixTransportPK(pk_list.sublist(pk.getName()), evals_list)

def initialConditionsList(xml):
    for pk in xml.sublist("PKs"):
        if pk.isElement("initial condition"):
            pk.sublist("initial condition").setName("initial conditions")


def update(xml):
    coupledFlowTransportMPC(xml)

    #enforceDtHistory(xml)
    timeStep(xml)
    tensorPerm(xml)
    initialConditionsList(xml)
    fixTransportPKs(xml)

    # this fixes chemistry
    fix_chemistry_ts_control.fixAll(xml)

    # moves lc types from ICs -> model parameters
    moveLC(xml)

    # viscosity relation type -> viscosity type
    viscosityRelationType(xml)

    # molar mass --> molar mass [kg m^-3]
    molarMass(xml)

    # tc-99 --> tracer
    removeTc99(xml)

    # verbose object controlled by global verbosity
    # NOTE: this is useful on tests, but a bad idea to use for user input files
    # removeVerboseObject(xml)
    
    

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
