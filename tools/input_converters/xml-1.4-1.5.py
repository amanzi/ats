#!/usr/bin/env python3
"""ATS input converter from version 1.4 to master branch"""

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
from convert_parameters_vg2bc import get_bc_param_from_vg


#
# helper functions
#
def rename_element(xml, old, new):
    try:
        child = asearch.child_by_name(xml, old)
    except aerrors.MissingXMLError:
        pass
    else:
        child.setName(new)

def retype_evaluator(xml, old_type, new_type):
    for e in asearch.find_path(xml, ["state", "evaluators"], no_skip=True):
        e_type = e.getElement("evaluator type")
        if e_type.getValue() == old_type:
            e_type.setValue(new_type)
        
def replace_string_in_name(xml, old, new):
    for match in xml.findall('.//'):
        if old in match.getName():
            match.setName(match.getName().replace(old, new))

def replace_string_in_value(xml, old, new):
    for match in xml.findall('.//'):
        if isinstance(match, parameter.Parameter) and match.getType() == 'string' and old in match.getValue():
            match.setValue(match.getValue().replace(old, new))


def _find_wrm_list(xml):
    """Finds the list, either in a flow PK or in a saturation evaluator"""
    try:
        return asearch.find_path(xml, ["PKs", "water retention evaluator"])
    except aerrors.MissingXMLError:
        fe_list = asearch.find_path(xml, ["state", "evaluators"], no_skip=True)
        try:
            return next(fe for fe in fe_list if fe.isElement("WRM parameters"))
        except StopIteration:
            try:
                return next(fe for fe in fe_list if fe.isElement("model parameters") and fe.getElement("model parameters").getValue() == "WRM parameters")
            except StopIteration:
                return None
    

def _has_surface_temp(sT):
    if sT.getElement("evaluator type").getValue() == "independent variable":
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
    """Deals with changes to surface temperatures."""
    # Many input files have a list where the surface temp was incorrectly set
    try:
        sT = asearch.find_path(xml, ["state", "evaluators", "surface-temperature"], no_skip=True)
    except aerrors.MissingXMLError:
        return

    if _has_surface_temp(sT):
        for func in asearch.findall_path(sT, ["function", "function", "function-composition"]):
            f1 = func.getElement("function1")
            f2 = func.getElement("function2")
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

    # remove the "air-snow temp diff" parameter in favor of creating a snow-temp field
    try:
        snowmelt = asearch.find_path(xml, ["state", "evaluators", "snow-melt"], no_skip=True)
    except aerrors.MissingXMLError:
        pass
    else:
        try:
            smr = snowmelt.getElement("air-snow temperature difference [C]")
        except aerrors.MissingXMLError:
            # already done?
            return

        try:
            snow_exp_temp = asearch.find_path(xml, ["state", "evaluators", "snow-expected_temperature"], no_skip=True)
        except aerrors.MissingXMLError:
            # make one
            eval_list = asearch.find_path(xml, ["state", "evaluators"], no_skip=True)
            snow_exp_temp = eval_list.sublist("snow-expected_temperature")
            snow_exp_temp.setParameter("evaluator type", 'string', "additive evaluator")
            snow_exp_temp.setParameter("dependencies", 'Array(string)', ["surface-air_temperature"])
            snow_exp_temp.setParameter("constant shift", 'double', -smr.getValue())
        else:
            pass

        try:
            snow_temp = asearch.find_path(xml, ["state", "evaluators", "snow-temperature"], no_skip=True)
        except aerrors.MissingXMLError:
            # make one
            eval_list = asearch.find_path(xml, ["state", "evaluators"], no_skip=True)
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

    # add a canopy-temperature field
    try:
        can_temp = asearch.find_path(xml, ["state", "evaluators", "canopy-temperature"], no_skip=True)
    except aerrors.MissingXMLError:
        # make one
        eval_list = asearch.find_path(xml, ["state", "evaluators"], no_skip=True)
        can_temp = eval_list.sublist("canopy-temperature")
        can_temp.setParameter("evaluator type", 'string', "additive evaluator")
        can_temp.setParameter("dependencies", 'Array(string)', ["surface-air_temperature"])
        can_temp.setParameter("constant shift", 'double', -1.5)

    # clean up references to these new variables
    for p in asearch.findall_path(xml, ["snow temperature key",]):
        if p.getValue() == "surface-temperature":
            p.setValue("snow-temperature")
    for p in asearch.findall_path(xml, ["canopy temperature key",]):
        if p.getValue() == "surface-temperature":
            p.setValue("canopy-temperature")

    # also snow-evap needs new snow tmep
    try:
        snow_evap = asearch.find_path(xml, ["state", "evaluators", "snow-evaporation"], no_skip=True)
    except aerrors.MissingXMLError:
        pass
    else:
        try:
            st = snow_evap.getElement("surface temperature key")
        except aerrors.MissingXMLError:
            snow_evap.setParameter("surface temperature key", "snow-temperature")
        else:
            if st.getValue() == "surface-temperature":
                st.setValue("snow-temperature")


def mafic_to_cap_pres(xml):
    """Converts mafic potential --> capillary pressure in land cover list"""
    try:
        lc_list = asearch.find_path(xml, ["state", "initial conditions", "land cover types"], no_skip=True)
    except aerrors.MissingXMLError:
        return
    for lc in lc_list:
        rename_element(lc, "mafic potential at fully closed stomata [Pa]",
                       "capillary pressure at fully closed stomata [Pa]")
        rename_element(lc, "mafic potential at fully open stomata [Pa]",
                       "capillary pressure at fully open stomata [Pa]")


def transpiration_relperm(xml):
    """Change the transpiration eval from the old default to the new default.

    Old was based on only TRF, new is based on both TRF and soil k_rel
    """
    retype_evaluator(xml, "transpiration distribution, rooting depth",
                     "transpiration distribution, relative permeability")

    try:
        asearch.find_path(xml, ["state", "evaluators"]).pop("plant_wilting_factor")
    except aerrors.MissingXMLError:
        pass


def add_soil_resistance(xml, rs_option="Sakagucki-Zeng"):
    """Soil resistance is now an evaluator, needed if using evap.

    Parameters
    ----------
    rs_option : string
      One of "Sakagucki-Zeng" or "Sellers"
    dessicated_zone_set : float, list, dict
      float: all soil types use the same dessicated_zone_thickness
      list: dessicated_zone_set = [] equal length with soil types
      dict: dessicated_zone_set = {key: soil_type, val: dessicated_zone_thickness}
    """
    eval_list = asearch.find_path(xml, ["state", "evaluators"], no_skip=True)
    try:
        seb_eval_list = next(e for e in eval_list if e.getElement("evaluator type").getValue() in 
                             ['surface energy balance, two components',
                              'surface energy balance, three components',
                              'evaporation downregulation, soil resistance'])
    except StopIteration:
        # no evap, so no soil resistance
        return

    # figure out the domain name
    domain = seb_eval_list.getName().split('-')[0] if '-' in seb_eval_list.getName() else ''
    if domain in ['', 'snow', 'surface']:
        sr_name = "surface-soil_resistance"
    else:
        if 'snow' in domain:
            domain = domain.replace('snow', 'surface')
        elif 'surface' in domain:
            pass
        else:
            domain = 'surface_'+domain
        sr_name = '-'.join([domain, 'soil_resistance'])

    try:
        rs_list = eval_list.getElement(sr_name)
    except aerrors.MissingXMLError:
        # no existing list, add one
        rs_list = eval_list.sublist(sr_name)
        rs_list.append(parameter.StringParameter("evaluator type", f'soil resistance, {rs_option}'))
        rs_list.append(parameter.StringParameter("model parameters", "WRM parameters"))

    if rs_option == 'Sakagucki-Zeng':
        # must also add dessicated zone thickness parameter
        wrm_param_list = asearch.find_path(xml, ["state", "model parameters", "WRM parameters"], no_skip=True)
        sub_list = asearch.children_by_tag(wrm_param_list, 'ParameterList')

        dzs = set(dz.getValue() for dz in asearch.findall_path(xml, ["state", "land cover types", "dessicated zone thickness [m]"]))
        if len(dzs) > 1:
            raise RuntimeError('Multiple dessicated zone thicknesses in land cover -- cannot map them to soil type automatically.')
        elif len(dzs) == 0:
            return
        else:
            dz_thick = dzs.pop()

        for s in sub_list:
            if not s.isElement("dessicated zone thickness [m]"):
                s.append(parameter.DoubleParameter("dessicated zone thickness [m]", dz_thick)) 



def add_wrm_to_model_parameters(xml):
    """Lifting WRM parameters from water retention in PK.

    Moves Richards PK-->water retention evaluator-->model parameters
    list to state->model parameters.
    """
    wrm_list = _find_wrm_list(xml)
    if wrm_list is None:
        return

    try:
        mpar = wrm_list.getElement("model parameters")
    except aerrors.MissingXMLError:
        wrm_list.append(parameter.StringParameter("model parameters", "WRM parameters"))
    else:
        mpar.setValue("WRM parameters")

    state_list = asearch.find_path(xml, ["state"])
    model_par_list = state_list.sublist("model parameters")
    if not model_par_list.isElement("WRM parameters"):
        try:
            wrm_param_list = asearch.remove_element(wrm_list, ["WRM parameters"], False, False)
        except aerrors.MissingXMLError:
            pass
        else:
            model_par_list.append(wrm_param_list)


def del_lc_params(xml):
    """Removes the now-dead land cover parameters for dessicated zone thickness and C&H b"""
    try:
        lc_list = asearch.find_path(xml, ["state", "initial conditions", "land cover types"], no_skip=True)
    except aerrors.MissingXMLError:
        return

    sub_lists = asearch.children_by_tag(lc_list, 'ParameterList')
    for s in sub_lists:
        try:
            s.pop("dessicated zone thickness [m]")
        except aerrors.MissingXMLError:
            pass

        try:
            s.pop("Clapp and Hornberger b [-]")
        except aerrors.MissingXMLError:
            pass
        
        

def add_rel_perm(xml, relp_option="relative permeability, water retention model"):
    """Ensures there is an evaluator for "relative_permeability"
    """
    eval_list = asearch.find_path(xml, ["state", "evaluators"], no_skip=True)

    try:
        perm_list = next(child for child in eval_list if
                         (child.getName() == "permeability" or 
                          child.getName().endswith("-permeability")))
    except StopIteration:
        return None # no perm, then no rel_perm at all!

    # get a prefix
    if perm_list.getName() == "permeability":
        relp_name = "relative_permeability"
    else:
        relp_name = perm_list.getName()[0:-len("permeability")]+"relative_permeability"

    try:
        relp_list = eval_list.getElement(relp_name)
    except aerrors.MissingXMLError:
        relp_list = eval_list.sublist(relp_name)
        relp_list.append(parameter.StringParameter("evaluator type", relp_option))
        if relp_option == "relative permeability, water retention model":
            relp_list.append(parameter.StringParameter("model parameters", "WRM parameters"))
        else:
            assert(relp_option == "relative permeability, freezing Brooks-Corey")
            relp_list.append(parameter.StringParameter("model parameters", "freezing rel perm parameters"))
            relp_list.append(parameter.DoubleParameter("omega [-]", 2))

    return relp_name


def move_wrm_relperm_params_to_relperm(xml, relp_name, par):
    try:
        relp_list = asearch.find_path(xml, ["state", "evaluators", relp_name])
    except aerrors.MissingXMLError:
        return

    if relp_list.isElement(par):
        return

    wrm_list = _find_wrm_list(xml)
    if wrm_list is None:
        return

    if wrm_list.isElement(par):
        relp_list.append(wrm_list.pop(par))




def move_wrm_to_state_evaluators(xml):
    """Last step -- moves the list out of Richards PK to state"""
    for pk in xml.getElement("PKs"):
        pk_type = pk.getElement("PK type").getValue()
        if pk_type in ["richards flow", "permafrost flow", "richards steady state"] and \
           pk.isElement("water retention evaluator"):
            wrm = pk.pop("water retention evaluator")
            wrm.setName("saturation_liquid")
            if pk_type == "permafrost flow":
                wrm.setParameter("evaluator type", "string", "water retention model with ice")
            else:
                wrm.setParameter("evaluator type", "string", "water retention model")

            wrm2 = copy.deepcopy(wrm)
            wrm2.setName("saturation_gas")
            evals_list = asearch.find_path(xml, ["state", "evaluators"], no_skip=True)
            evals_list.append(wrm)
            evals_list.append(wrm2)

    
def add_frz_relp_to_model_parameters(xml):
    def copy_element(a_soil_list, var_name, var_type):
        if var_type == "string":
            return parameter.StringParameter(var_name, a_soil_list.getElement(var_name).getValue())
        elif var_type == "double":
            return parameter.DoubleParameter(var_name, a_soil_list.getElement(var_name).getValue())
        else:
            pass

    model_par_list = asearch.find_path(xml, ["state", "model parameters"], no_skip=True)
    frz_relp_par_list = model_par_list.sublist("freezing rel perm parameters")
    try:
        soil_types = asearch.children_by_tag(asearch.remove_element(xml, 
                ["water retention evaluator", "WRM parameters"], 
                False, False), 'ParameterList')
    except aerrors.MissingXMLError:
        soil_types = asearch.children_by_tag(model_par_list.getElement("WRM parameters"), 
                                             'ParameterList')
    for s in soil_types:
        slist = frz_relp_par_list.sublist(s.getName())
        for var_name, var_type in zip(["region", "residual saturation [-]", 
            "smoothing interval width [saturation]"], ["string", "double", "double"]):
            slist.append(copy_element(s, var_name, var_type))

        slist.append(parameter.StringParameter("WRM Type", "Brooks-Corey"))
        alpha = s.getElement("van Genuchten alpha [Pa^-1]").getValue()
        try:
            n = s.getElement("van Genuchten n [-]").getValue()
        except:
            m = s.getElement("van Genuchten m [-]").getValue()
            n = 1 / (1 - m)
        bc_satp, bc_lambda_recip = get_bc_param_from_vg(alpha, n)
        slist.append(parameter.DoubleParameter("Brooks-Corey saturated matric suction [Pa]", bc_satp))
        slist.append(parameter.DoubleParameter("Brooks-Corey lambda [-]", 1 / bc_lambda_recip))


def lowercase_wrmtype(xml):
    replace_string_in_name(xml, "WRM type", "wrm type")
    replace_string_in_name(xml, "WRM Type", "wrm type")
    replace_string_in_name(xml, "permafrost WRM type", "permafrost wrm type")
    replace_string_in_name(xml, "permafrost WRM Type", "permafrost wrm type")

    
def brooks_corey(xml):
    replace_string_in_name(xml, "Brooks Corey", "Brooks-Corey")
    replace_string_in_value(xml, "Brooks Corey", "Brooks-Corey")

            
def update(xml, transpiration_distribution=False, frozen_krel=False,
           soil_res='Sakagucki-Zeng'):
    """generic update calls all needed things"""
    mafic_to_cap_pres(xml)
    replace_string_in_name(xml, "rooting_depth_fraction", "root_fraction")
    retype_evaluator(xml, "rooting depth fraction", "root fraction")
    surface_temp(xml)

    retype_evaluator(xml, "transpiration distribution via rooting depth", "transpiration distribution, rooting depth")
    retype_evaluator(xml, "transpiration distribution via relative permeability", "transpiration distribution, relative permeability")
    retype_evaluator(xml, "evaporation downregulation via soil resistance", "evaporation downregulation, soil resistance")

    if transpiration_distribution:
        transpiration_relperm(xml)

    add_wrm_to_model_parameters(xml)
    add_soil_resistance(xml, soil_res)
    del_lc_params(xml)

    if frozen_krel:
        relp_name = add_rel_perm(xml, "relative permeability, freezing Brooks-Corey")
        add_frz_relp_to_model_parameters(xml)
    else:
        relp_name = add_rel_perm(xml, "relative permeability, water retention model")

    if relp_name is not None:
        for par in ["minimum rel perm cutoff", "use surface rel perm"]:
            move_wrm_relperm_params_to_relperm(xml, relp_name, par)

    move_wrm_to_state_evaluators(xml)
    retype_evaluator(xml, "WRM rel perm", "relative permeability, water retention model")
    retype_evaluator(xml, "relative permeability, van Genuchten", "relative permeability, water retention model")
    retype_evaluator(xml, "WRM", "water retention model")
    retype_evaluator(xml, "permafrost WRM", "water retention model with ice")
    lowercase_wrmtype(xml)
    brooks_corey(xml)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("infile", help="input filename")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-i", "--inplace", action="store_true", help="fix file in place")
    group.add_argument("-o", "--outfile", help="output filename")

    parser.add_argument("--transpiration-distribution", action="store_true",
                        help="Update to the relperm based transpiration distribution")
    parser.add_argument("--frozen-krel", action="store_true",
                        help="Update relperm to the new Arctic variant")
    parser.add_argument("--soil-resistance", default='Sakagucki-Zeng', choices=['Sakagucki-Zeng', 'Sellers'],
                        help="Soil resistance model name")
    
    args = parser.parse_args()

    # check for orig file
    print("Converting file: %s"%args.infile)
    xml = aio.fromFile(args.infile, True)
    update(xml, args.transpiration_distribution, args.frozen_krel,
           args.soil_resistance)
    if args.inplace:
        aio.toFile(xml, args.infile)
    else:
        aio.toFile(xml, args.outfile)
    sys.exit(0)
    

    
