#!/usr/bin/env python3
"""ATS input converter from 1.0 to 1.1"""

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
    fe_list = asearch.find_path(xml, ["state","field evaluators"])
    for eval in fe_list:
        ename = eval.getName()
        if ename.endswith('manning_coefficient') and not ename.startswith('snow'):
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
                            assert len(asearch.findall_name(reg, 'components')) == 0
                            reg.append(parameter.ArrayStringParameter('components', ['cell', 'boundary_face']))
                            fixed = True

                    if not fixed:
                        try:
                            comp = reg.getElement('components')
                        except aerrors.MissingXMLError:
                            pass
                        else:
                            reg.pop('components')
                            assert len(asearch.findall_name(reg, 'components')) == 0
                            reg.append(parameter.ArrayStringParameter('components', ['cell', 'boundary_face']))
                            fixed = True

                    if not fixed:
                        print(reg.__str__().decode('utf-8'))
                        raise aerrors.MissingXMLError('Missing "component" or "components"')


def fixSnow(xml):
    # remove "include dt factor"
    # rename "dt factor" --> "dt factor [s]"
    try:
        snow_cond = asearch.find_path(xml, ["state", "field evaluators", "snow-conductivity"])
    except aerrors.MissingXMLError:
        pass
    else:
        try:
            snow_cond.pop("include dt factor")
        except aerrors.MissingXMLError:
            pass

        try:
            dtf = asearch.find_name(snow_cond, "dt factor")
        except aerrors.MissingXMLError:
            pass
        else:
            dtf.setName("dt factor [s]")

        # rename "include density factor" --> "include density"
        try:
            inc_dens = asearch.find_name(snow_cond, "include density factor")
        except aerrors.MissingXMLError:
            pass
        else:
            inc_dens.setName("include density")

        # add "swe density factor [-]" (default is 10)
        if len(asearch.findall_name(snow_cond, "swe density factor [-]")) == 0:
            snow_cond.append(parameter.DoubleParameter("swe density factor [-]", 10.0))

            
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
            

def _get_prefixed_name(pk, name):
    try:
        domain = asearch.find_name(pk, "domain name").getValue()
    except aerrors.MissingXMLError:
        domain = "domain"

    if domain == "domain":
        default_val = name
    else:
        default_val = domain+"-"+name
    return default_val    

                
def fixMassSource(xml):
    # find all pks
    pks = asearch.find_name(xml, "PKs")
    for pk in pks:
        pk_type = asearch.find_name(pk, "PK type")
        if pk_type.getValue() in ["richards",
                                  "overland flow, pressure basis",
                                  "overland flow",
                                  "overland flow with ice",
                                  "snow distribution",
                                  "richards steady state",
                                  "permafrost flow"]:
            # has source?
            try:
                has_source = asearch.find_name(pk, "source term")
            except aerrors.MissingXMLError:
                continue

            if has_source.getValue():
                # check for the default key
                try:
                    source_name = asearch.find_name(pk, "source key")
                except aerrors.MissingXMLError:
                    try:
                        source_name_suffix = asearch.find_name(pk, "source key suffix")
                    except aerrors.MissingXMLError:
                        # if we are using the default, check to see if
                        # the old default is in the list of evaluators
                        default_val = _get_prefixed_name(pk, "mass_source")
                        
                        try:
                            source_eval = asearch.find_path(xml, ["state","field evaluators",default_val])
                        except aerrors.MissingXMLError:
                            # likely this is a Arctic run, and the
                            # source name is now changed to the new
                            # default, do nothing
                            pass
                        else:
                            # we are using the default, change it to the new default
                            prefix_name = _get_prefixed_name(pk, "water_source")
                            source_eval.setName(prefix_name)
                    else:
                        if source_name_suffix.getValue() == "mass_source":
                            # if we found source key suffix, and it is mass_source...
                            default_val = _get_prefixed_name(pk, "mass_source")
                            try:
                                source_eval = asearch.find_path(xml, ["state","field evaluators",default_val])
                            except aerrors.MissingXMLError:
                                # likely this is an Arctic run, and
                                # the source name is now changed to
                                # the new default.  Change the suffix
                                source_name_suffix.setValue("water_source")
                            else:
                                # change both the suffix and the evaluator
                                source_name_suffix.setValue("water_source")
                                prefix_name = _get_prefixed_name(pk, "water_source")
                                source_eval.setName(_get_prefixed_name(pk, "water_source"))

                        
                else:
                    if source_name.getValue().endswith("mass_source"):
                        # we found a source key, and its suffix is mass source...
                        try:
                            source_eval = asearch.find_path(xml, ["state","field evaluators",source_name.getValue()])
                        except aerrors.MissingXMLError:
                            # likely Arctic
                            if '-' in source_name.getValue():
                                source_domain = source_name.getValue().split('-')[0]
                                source_name.setValue(source_domain+'-water_source')
                            else:
                                source_name.setValue('water_source')
                        else:
                            source_domain = source_name.getValue().split('-')[0]
                            source_name.setValue(source_domain+'-water_source')
                            source_eval.setName(source_domain+'-water_source')

                    
    # changes "mass source in meters" to "water source".  Note this
    # can be in transport, not only in flow.
    for ws_in_meters in asearch.findall_name(xml, "mass source in meters"):
        ws_in_meters.setName("water source in meters")

    # changes "mass source units" to "water source units" in total energy sources
    fe_list = asearch.find_path(xml, ["state","field evaluators"])
    for eval in fe_list:
        if asearch.find_name(eval, "field evaluator type").getValue() == "advected energy source":
            try:
                units = asearch.find_name(eval, "mass source units")
            except aerrors.MissingXMLError:
                pass
            else:
                units.setName("water source units")

        
def update(xml):
    vanGenuchtenParams(xml)
    thermalConductivityParams(xml)
    commonSupressOptions(xml)
    
    checkManning(xml)
    # NOTE: these will get added when subgrid pull request is done
    fixSnow(xml)
    #fixSubgrid(xml)
    mergePreconditionerLinearSolver(xml)
    fixOverlandConductivity(xml)
    updateObservations(xml)
    fixMassSource(xml)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Fix a number of changes from ATS input spec 1.0 to 1.1")
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
    

    
