#!/usr/bin/env python

"""ATS input converter from 0.86 to 0.87"""

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

import flatten_pks

def fixEvaluator(xml, name, newname):
    try:
        pd = asearch.find_path(xml, ["state","field evaluators",name])
    except aerrors.MissingXMLError:
        pass
    else:
        pd.setName(newname)
    

def compressibility(xml):
    for c in asearch.findall_name(xml, "pore compressibility"):
        c.set("name", "pore compressibility [Pa^-1]")

def diffusion(xml):
    for diff in asearch.findall_name(xml, "Diffusion"):
        diff.set("name", "diffusion")
    for diff in asearch.findall_name(xml, "Diffusion PC"):
        diff.set("name", "diffusion preconditioner")

def water_energy(xml):
    try:
        water = asearch.find_path(xml, ["state","field evaluators","water_content"])
    except aerrors.MissingXMLError:
        pass
    else:
        ftype = water.getElement("field evaluator type")
        if ftype.value == "permafrost water content":
            if water.isElement("include water vapor") and water.getElement("include water vapor").value:
                ftype.setValue("three phase water content")
            else:
                ftype.setValue("liquid+ice water content")
                    
        elif ftype.value == "richards water content":
            if water.isElement("include water vapor") and water.getElement("include water vapor").value:
                ftype.setValue("liquid+gas water content")
            else:
                ftype.setValue("richards water content")

    try:
        energy = asearch.find_path(xml, ["state","field evaluators","energy"])
    except aerrors.MissingXMLError:
        pass
    else:
        etype = energy.getElement("field evaluator type")
        if ftype.value == "three phase water content":
            etype.setValue("three phase energy")
        elif ftype.value == "liquid+ice water content":
            etype.setValue("liquid+ice energy")
        elif ftype.value == "liquid+gas water content":
            etype.setValue("liquid+gas energy")
        elif ftype.value == "richards water content":
            etype.setValue("richards energy")
            
def adds_source_units(xml):
    fevals = asearch.find_path(xml, ["state","field evaluators"])
    try:
        te = fevals.getElement("total_energy_source")
    except aerrors.MissingXMLError:
        pass
    else:
        if not te.isElement("mass source units"):
            te.setParameter("mass source units", "string", "mol m^-2 s^-1")

    try:
        tes = fevals.getElement("surface-total_energy_source")
    except aerrors.MissingXMLError:
        pass
    else:
        if not tes.isElement("mass source units"):
            tes.setParameter("mass source units", "string", "m s^-1")


def seepage_face_bcs(xml):
    for bclist in asearch.findall_name(xml, "boundary conditions"):
        if bclist.isElement("seepage face"):
            agroup = bclist.sublist("seepage face")[0]
            if agroup.isElement("boundary pressure"):
                bclist.sublist("seepage face").set("name", "seepage face pressure")
            elif agroup.isElement("boundary head"):
                bclist.sublist("seepage face").set("name", "seepage face head")
            else:
                raise RuntimeError("unknown seepage condition")

def primary_variable(xml):
    for pv in asearch.findall_name(xml, "primary variable"):
        pv.setName("primary variable key")
    for pv in asearch.findall_name(xml, "primary variable key"):
        if pv.getValue() == "ponded_depth":
            pv.setValue("surface-ponded_depth")
        if pv.getValue() == "snow_depth":
            pv.setValue("surface-snow_depth")
            
def mesh_list(xml):
    try:
        xml.pop("Native Unstructured Input")
    except aerrors.MissingXMLError:
        pass

    try:
        xml.pop("grid_option")
    except aerrors.MissingXMLError:
        pass

    # move domain mesh parameters to a sublist
    mesh = asearch.child_by_name(xml, "mesh")
    domain = mesh.sublist("domain")
    to_pop = []
    for el in mesh:
        if el.get("name") == "surface mesh":
            el.setName("surface")
        elif el.get("name") == "column meshes":
            el.setName("column")
        elif el.get("name") == "column surface meshes":
            el.setName("column surface")
        elif el.get("name") in ["surface", "column", "column surface", "domain", "subgrid"]:
            pass
        elif el.get("name") in ["framework"]:
            to_pop.append(el)
        else:
            domain.append(el)

    if domain.isElement("framework"):
        domain.pop("framework")
            
    for el in domain:
        if mesh.isElement(el.get("name")):
            mesh.pop(el.get("name"))
    for el in to_pop:
        mesh.pop(el.get("name"))        

    # move surface mesh parameters to a sublist
    if mesh.isElement("surface"):
        surf_list = mesh.sublist("surface")
        surf_p_list = surf_list.sublist("surface")
        if surf_list.isElement("surface sideset name"):
            surf_p_list.append(surf_list.pop("surface sideset name"))
        if surf_list.isElement("surface sideset names"):
            surf_p_list.append(surf_list.pop("surface sideset names"))
        
    # make sure all left are mesh sublists, add a mesh type parameter
    valid_types = ["read mesh file", "generate mesh", "logical mesh",
                   "aliased", "surface", "column", "column surface", "subgrid"]    
    for el in mesh:
        if not el.isElement("mesh type"):
            print(el)
            found = False
            for valid_type in valid_types:
                if el.isElement(valid_type):
                    print("setting type: ", valid_type)
                    el.setParameter("mesh type", "string", valid_type)
                    asearch.child_by_name(el, valid_type).setName(valid_type+" parameters")
                    found = True
                    break
            assert(found)
        
    return


def vis(xml):
    if xml.isElement("visualization") and not asearch.child_by_name(xml, "visualization").isElement("domain"):
        vis_domain = xml.pop("visualization")
        vis_domain.setName("domain")
    
        vis_list = xml.sublist("visualization")
        vis_list.append(vis_domain)
        if xml.isElement("visualization surface"):
            vis_surf = xml.pop("visualization surface")
            vis_surf.setName("surface")
            vis_list.append(vis_surf)

        if xml.isElement("visualization columns"):
            vis_col = xml.pop("visualization columns")
            vis_col.setName("column_*")
            vis_list.append(vis_col)


def snow_depth(xml):
    for seb_pk in asearch.gen_by_path(xml, ["PKs","SEB"]):
        if seb_pk.isElement("primary variable key") and \
             asearch.child_by_name(seb_pk,"primary variable key").getValue() == "surface-snow_depth" and \
             seb_pk.isElement("conserved quantity key") and \
             asearch.child_by_name(seb_pk,"conserved quantity key").getValue() == "snow_depth":
            asearch.child_by_name(seb_pk,"conserved quantity key").set("value", "surface-snow_depth")

def snow_distribution(xml):
    for snow_dist_pk in asearch.findall_path(xml, ["PKs","snow distribution"]):
        snow_dist_pk.append(parameter.DoubleParameter("distribution time", 86400.0))
        if snow_dist_pk.isElement("primary variable key") and \
             asearch.child_by_name(snow_dist_pk,"primary variable key").getValue() == "surface-precipitation_snow" and \
             snow_dist_pk.isElement("conserved quantity key") and \
             asearch.child_by_name(snow_dist_pk,"conserved quantity key").getValue() == "precipitation-snow":
            asearch.child_by_name(snow_dist_pk,"conserved quantity key").set("value", "surface-precipitation-snow")
        
    for ssk in asearch.findall_path(xml, ["state","field evaluators","surface-snow_skin_potential"]):
        if not ssk.isElement("dt factor"):
            ssk.append(parameter.DoubleParameter("dt factor", 86400.0))
        else:
            asearch.child_by_name(ssk, "dt factor").set("value","86400.0")

    for ssc in asearch.findall_path(xml, ["state", "field evaluators","surface-snow_conductivity"]):
        if not ssc.isElement("include dt factor"):
            ssc.append(parameter.BoolParameter("include dt factor", True))
        else:
            asearch.child_by_name(ssc, "include dt factor").set("value","true")

        if not ssc.isElement("dt factor"):
            ssc.append(parameter.DoubleParameter("dt factor", 86400.0))
        else:
            asearch.child_by_name(ssc, "dt factor").set("value","86400.0")


def seb(xml):
    pk_list = asearch.findall_name(xml, "PKs")
    assert(len(pk_list) == 1)
    pk_list = pk_list[0]

    seb_pk = None
    flow_sub_pk = None
    flow_surf_pk = None
    energy_sub_pk = None
    energy_surf_pk = None
    # make sure we can find all of them!
    for pk in pk_list:
        if asearch.child_by_name(pk,"PK type").getValue() == "surface balance implicit":
            if seb_pk is not None:
                raise RuntimeError("Cannot deal with SEB changes!")
            seb_pk = pk
        elif asearch.child_by_name(pk,"PK type").getValue() == "permafrost flow":
            if flow_sub_pk is not None:
                raise RuntimeError("Cannot deal with SEB changes!")
            flow_sub_pk = pk
        elif asearch.child_by_name(pk,"PK type").getValue() == "overland flow with ice":
            if flow_surf_pk is not None:
                raise RuntimeError("Cannot deal with SEB changes!")
            flow_surf_pk = pk
        elif asearch.child_by_name(pk,"PK type").getValue() == "three-phase energy":
            if energy_sub_pk is not None:
                raise RuntimeError("Cannot deal with SEB changes!")
            energy_sub_pk = pk
        elif asearch.child_by_name(pk,"PK type").getValue() == "surface energy":
            if energy_surf_pk is not None:
                raise RuntimeError("Cannot deal with SEB changes!")
            energy_surf_pk = pk

    if seb_pk is None or flow_sub_pk is None or flow_surf_pk is None or energy_sub_pk is None or energy_surf_pk is None:
        return

    # check the source terms for all
    def set_source_term(pk):
        if not pk.isElement("source term"):
            pk.append(parameter.BoolParameter("source term", True))
        else:
            asearch.child_by_name(pk, "source term").set("value","true")
    set_source_term(flow_sub_pk)
    set_source_term(flow_surf_pk)
    set_source_term(energy_sub_pk)
    set_source_term(energy_surf_pk)

    if not flow_sub_pk.isElement("mass source key"):
        flow_sub_pk.append(parameter.StringParameter("mass source key", "mass_source"))
    if not flow_surf_pk.isElement("source key"):
        flow_surf_pk.append(parameter.StringParameter("source key", "surface-mass_source"))
    if not flow_surf_pk.isElement("mass source in meters"):
        flow_surf_pk.append(parameter.BoolParameter("mass source in meters", True))
    if not energy_sub_pk.isElement("energy source"):
        energy_sub_pk.append(parameter.StringParameter("energy source", "total_energy_source"))
    if not energy_surf_pk.isElement("energy source"):
        energy_surf_pk.append(parameter.StringParameter("energy source", "surface-total_energy_source"))

    eval_list = asearch.find_path(xml, ["state","field evaluators"])
    try:
        eval_list.pop("surface-total_energy_source")
    except aerrors.MissingXMLError:
        pass
    try:
        eval_list.pop("surface-mass_source_enthalpy")
    except aerrors.MissingXMLError:
        pass
    try:
        eval_list.pop("surface-source_internal_energy")
    except aerrors.MissingXMLError:
        pass

    molar_dens = asearch.child_by_name(eval_list, "surface-source_molar_density")
    if molar_dens.isElement("temperature key"):
        asearch.child_by_name(molar_dens, "temperature key").set("value", "surface-temperature")
    else:
        molar_dens.append(parameter.StringParameter("temperature key", "surface-temperature"))
    
            
            
def update(xml):
    flatten_pks.flatten_pks(xml)
    fixEvaluator(xml, "ponded_depth", "surface-ponded_depth")
    fixEvaluator(xml, "ponded_depth_bar", "surface-ponded_depth_bar")
    fixEvaluator(xml, "manning_coefficient", "surface-manning_coefficient")
    fixEvaluator(xml, "unfrozen_fraction", "surface-unfrozen_fraction")
    fixEvaluator(xml, "unfrozen_effective_depth", "surface-unfrozen_effective_depth")
    fixEvaluator(xml, "incoming_shortwave_radiation", "surface-incoming_shortwave_radiation")
    fixEvaluator(xml, "co2_concentration", "surface-co2_concentration")
    fixEvaluator(xml, "precipitation_snow", "surface-precipitation_snow")
    fixEvaluator(xml, "precipitation_rain", "surface-precipitation_rain")
    fixEvaluator(xml, "relative_humidity", "surface-relative_humidity")
    fixEvaluator(xml, "wind_speed", "surface-wind_speed")
    fixEvaluator(xml, "incoming_longwave_radiation", "surface-incoming_longwave_radiation")
    fixEvaluator(xml, "air_temperature", "surface-air_temperature")
    compressibility(xml)
    diffusion(xml)
    water_energy(xml)
    adds_source_units(xml)
    seepage_face_bcs(xml)
    primary_variable(xml)
    #mesh_list(xml)
    vis(xml)
    snow_depth(xml)
    snow_distribution(xml)
    seb(xml)

    import verbose_object
    verbose_object.fixVerboseObject(xml)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Fix a number of changes from ATS input spec 0.86 to 0.9x")
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
    

    
