#!/usr/bin/env python3
"""ATS input converter from 0.88 to master"""

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
from amanzi_xml.common import parameter

def change_name(xml, old, new, allow_multiple=False):
    try:
        asearch.change_name(xml, old, new, allow_multiple)
    except aerrors.MissingXMLError:
        pass


def vanGenuchtenParams(xml):
    """Can turn off derivative of source terms"""
    for wrm in asearch.gen_by_path(xml, ["WRM parameters"]):
        for region in wrm:
            change_name(region, "van Genuchten alpha", "van Genuchten alpha [Pa^-1]")
            change_name(region, "van Genuchten m", "van Genuchten m [-]")
            change_name(region, "van Genuchten n", "van Genuchten n [-]")
            change_name(region, "van Genuchten residual saturation", "residual saturation [-]")
            change_name(region, "residual saturation", "residual saturation [-]")

def thermalConductivityParams(xml):
    change_name(xml, "thermal conductivity of soil [W/(m-K)]", "thermal conductivity of soil [W m^-1 K^-1]", allow_multiple=True)
    change_name(xml, "thermal conductivity of ice [W/(m-K)]", "thermal conductivity of ice [W m^-1 K^-1]", allow_multiple=True)
    change_name(xml, "thermal conductivity of liquid [W/(m-K)]", "thermal conductivity of liquid [W m^-1 K^-1]", allow_multiple=True)
    change_name(xml, "thermal conductivity of water [W/(m-K)]", "thermal conductivity of water [W m^-1 K^-1]", allow_multiple=True)
    change_name(xml, "thermal conductivity of gas [W/(m-K)]", "thermal conductivity of gas [W m^-1 K^-1]", allow_multiple=True)
    change_name(xml, "thermal conductivity of frozen zone [W/(m-K)]", "thermal conductivity of frozen zone [W m^-1 K^-1]", allow_multiple=True)
    change_name(xml, "thermal conductivity of unfrozen zone [W/(m-K)]", "thermal conductivity of unfrozen zone [W m^-1 K^-1]", allow_multiple=True)
    change_name(xml, "thermal conductivity of mushy zone [W/(m-K)]", "thermal conductivity of mushy zone [W m^-1 K^-1]", allow_multiple=True)
    change_name(xml, "thermal conductivity, dry [W/(m-K)]", "thermal conductivity, dry [W m^-1 K^-1]", allow_multiple=True)
    change_name(xml, "thermal conductivity, wet [W/(m-K)]", "thermal conductivity, wet [W m^-1 K^-1]", allow_multiple=True)
    change_name(xml, "thermal conductivity, saturated (unfrozen) [W/(m-K)]",
                "thermal conductivity, saturated (unfrozen) [W m^-1 K^-1]", allow_multiple=True)
    change_name(xml, "quadratic u_0 [J/kg]", "latent heat [J kg^-1]", True)
    change_name(xml, "quadratic u_0 [J/mol]", "latent heat [J mol^-1]", True)
    change_name(xml, "quadratic a [J/kg-K]", "heat capacity [J kg^-1 K^-1]", True)
    change_name(xml, "quadratic a [J/mol-K]", "heat capacity [J mol^-1 K^-1]", True)
    change_name(xml, "quadratic b [J/kg-K^2]", "quadratic b [J kg^-1 K^-2]", True)
    change_name(xml, "quadratic b [J/mol-K^2]", "quadratic b [J mol^-1 K^-2]", True)
    change_name(xml, "latent heat [J/kg]", "latent heat [J kg^-1]", True)
    change_name(xml, "latent heat [J/mol]", "latent heat [J mol^-1]", True)
    change_name(xml, "heat capacity [J/kg-K]", "heat capacity [J kg^-1 K^-1]", True)
    change_name(xml, "heat capacity [J/mol-K]", "heat capacity [J mol^-1 K^-1]", True)
    change_name(xml, "Reference temperature [K]", "reference temperature [K]", True)
    change_name(xml, "heat capacity of air [J/(mol-K)]", "heat capacity [J mol^-1 K^-1]", True)
    change_name(xml, "heat of vaporization of water [J/mol]", "latent heat [J mol^-1]", True)
    change_name(xml, "capillary pressure [J/mol]", "latent heat [J mol^-1]", True)
    change_name(xml, "interfacial tension ice-water", "interfacial tension ice-water [mN m^-1]", True) # note flip in order of phases!
    change_name(xml, "interfacial tension air-water", "interfacial tension air-water [mN m^-1]", True)
    change_name(xml, "heat of fusion reference temperature [K]", "reference temperature [K]", True)
    change_name(xml, "heat of fusion of water [J/mol]", "latent heat [J mol^-1]", True)
    change_name(xml, "heat of fusion of water [J/kg]", "latent heat [J kg^-1]", True)
    change_name(xml, "smoothing length [K]", "smoothing width [K]", True)

def commonSupressOptions(xml):
    change_name(xml, "supress Jacobian terms: div hq / dp,T", "supress Jacobian terms: d div hq / dp,T")
    
            
def checkManning(xml):
    fe_list = asearch.gen_by_path(xml, ["state","field evaluators"])
    for eval in fe_list:
        if eval.get("name").endswith('manning_coefficient'):
            eval_type = eval.getElement('field evaluator type')
            if eval_type.getValue() == 'independent variable':
                func_reg = eval.getElement("function")
                for reg in func_reg:
                    comp_entries = asearch.findall_name(reg, ['components'])
                    if len(comp_entries) > 1:
                        # previous iterations of this script were broken...
                        for entry in comp_entries[1:]:
                            reg.remove(entry)
                        

                    fixed = False
                    if not fixed:
                        try:
                            comp = reg.getElement('component')
                        except aerrors.MissingXMLError:
                            pass
                        else:
                            reg.pop('component')
                            assert not any(el.getName() == 'component' for el in reg)
                            reg.append(parameter.ArrayStringParameter('components', ['cell', 'boundary_face']))
                            fixed = True

                    if not fixed:
                        try:
                            comp = reg.getElement('components')
                        except aerrors.MissingXMLError:
                            pass
                        else:
                            reg.pop('components')
                            assert not any(el.getName() == 'components' for el in reg)
                            reg.append(parameter.ArrayStringParameter('components', ['cell', 'boundary_face']))
                            fixed = True

                    if not fixed:
                        print(reg.__str__().decode('utf-8'))
                        raise aerrors.MissingXMLError('Missing "component" or "components"')


def fixSnow(xml):
    # todo: remove "include dt factor"
    # todo: rename "dt factor" --> "dt factor [s]"
    # todo: add "swe density factor [-]" (default is 10)
    # todo: rename "include density factor" --> "include density"
    raise NotImplementedError("fix snow")

def fixSubgrid(xml):
    # todo: ponded_depth_minus_depression_depth --> mobile_depth
    raise NotImplementedError("fix subgrid")

def mergePreconditionerLinearSolver(xml):
    pc_list = xml.getElement("PKs")
    for pk in pc_list:
        if pk.isElement("preconditioner"):
            pc = pk.getElement("preconditioner")
            pc.setName("inverse")
            change_name(pc, "preconditioner type", "preconditioning method")
            change_name(pc, "preconditioner method", "preconditioning method")
            if pk.isElement("linear solver"):
                ls = pk.getElement("linear solver")
                pc.extend(list(pk.getElement("linear solver")))
                pk.remove(ls)
                
        elif pk.isElement("linear solver"):
            ls = pk.getElement("linear solver")
            ls.setName("inverse")

        elif pk.isElement("inverse"):
            # already done, but check that preconditioning method is set
            pc = pk.getElement("inverse")
            change_name(pc, "preconditioner type", "preconditioning method")
            change_name(pc, "preconditioner method", "preconditioning method") 

def fixOverlandConductivity(xml):
    pm = asearch.parent_map(xml)
    fe = asearch.find_path(xml, ["state", "field evaluators"])
    try:
        elev = asearch.find_name(xml, "elevation evaluator")
    except aerrors.MissingXMLError:
        pass
    else:
        pk = pm[elev]
        pk.remove(elev)
        domain_name = pk.getElement("domain name").getValue()
        elev.setName(domain_name+"-elevation")
        if elev.isElement("elevation function"):
            elev.setParameter("field evaluator type", "string", "standalone elevation")
            for comp in asearch.findall_name(elev, "components"):
                comp.setValue(["cell", "face"])
        else:
            elev.setParameter("field evaluator type", "string", "meshed elevation")
        fe.append(elev)

    try:
        oc = asearch.find_name(xml, "overland conductivity evaluator")
    except aerrors.MissingXMLError:
        pass
    else:
        pk = pm[oc]
        pk.remove(oc)
        domain_name = pk.getElement("domain name").getValue()
        oc.setName(domain_name+"-overland_conductivity")
        oc.setParameter("field evaluator type", "string", "overland conductivity")
        fe.append(oc)

def updateObservations(xml):
    try:
        obs = asearch.findall_path(xml, ["observations", "functional"])
    except aerrors.MissingXMLError:
        pass
    else:
        for func in obs:
            func_val = func.getValue()
            if func_val.startswith("observation data: "):
                func.setValue(func_val[len("observation data: "):])
            

def update(xml):
    vanGenuchtenParams(xml)
    thermalConductivityParams(xml)
    commonSupressOptions(xml)
    
    checkManning(xml)
    # NOTE: these will get added when subgrid pull request is done
    #fixSnow(xml)
    #fixSubgrid(xml)
    mergePreconditionerLinearSolver(xml)
    fixOverlandConductivity(xml)
    updateObservations(xml)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Fix a number of changes from ATS input spec 0.88 to master")
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
    

    
