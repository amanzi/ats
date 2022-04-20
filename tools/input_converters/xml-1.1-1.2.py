#!/usr/bin/env python3
"""ATS input converter from 1.1 to release 1.2"""

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

def alias_mesh_target(xml):
    mesh_list = asearch.child_by_name(xml, "mesh")
    for mesh in mesh_list:
        mtype = asearch.child_by_name(mesh, "mesh type")
        if mtype.getValue() == "aliased":
            try:
                alias = asearch.find_path(mesh, ["aliased parameters", "alias"])
            except aerrors.MissingXMLError:
                pass
            else:
                alias.setName("target")


def seb_twocomponent(xml, eval_name):
    """This changes the old, 'Arctic' SEB model into a test-reproducing, same answer input file."""
    # area fractions
    try:
        frac_area = asearch.find_path(xml, ["state","field evaluators", "snow-fractional_areas"])
    except aerrors.MissingXMLError:
        pass
    else:
        frac_area.setName("surface-fractional_areas")

    for name in ['surface-fractional_areas', 'surface_column:*-fractional_areas']:
        try:
            frac_area = asearch.find_path(xml, ["state","field evaluators", name])
        except aerrors.MissingXMLError:
            pass
        else:
            domain = frac_area.getName().split('-')[0]
            frac_area.setName(domain+"-area_fractions")
            frac_area_type = asearch.child_by_name(frac_area, "field evaluator type")
            if frac_area_type.getValue() == "surface balance area fractions":
                frac_area_type.setValue("area fractions, two components")

    # water source
    try:
        water_source = asearch.find_path(xml, ['state', 'field evaluators', eval_name])
    except aerrors.MissingXMLError:
        pass
    else:
        water_source_type = asearch.child_by_name(water_source, "field evaluator type")
        if water_source_type.getValue() == "surface balance":
            water_source_type.setValue("surface energy balance, two components")

            try:
                asearch.child_by_name(water_source, "use model from ATS 1.1")
            except aerrors.MissingXMLError:
                water_source.append(parameter.BoolParameter("use model from ATS 1.1", True))

    # subgrid albedos
    for name in ['surface-subgrid_albedos', 'surface_column:*-subgrid_albedos']:
        try:
            albedo = asearch.find_path(xml, ['state', 'field evaluators', name])
        except aerrors.MissingXMLError:
            pass
        else:
            dname = name.split('-')[0]
            albedo.setName(dname+'-albedos')
            albedo_type = asearch.child_by_name(albedo, "field evaluator type")
            if albedo_type.getValue() == "albedo":
                albedo_type.setValue("subgrid albedos, two components")



def seb_threecomponent(xml, eval_name):
    """This changes the old, 'Arctic' SEB model into more robust variant that crashes less often."""
    # area fractions
    try:
        frac_area = asearch.find_path(xml, ["state","field evaluators", "snow-fractional_areas"])
    except aerrors.MissingXMLError:
        pass
    else:
        frac_area.setName("surface-fractional_areas")

    for name in ['surface-fractional_areas', 'surface_column:*-fractional_areas']:
        try:
            frac_area = asearch.find_path(xml, ["state","field evaluators", name])
        except aerrors.MissingXMLError:
            pass
        else:
            domain = frac_area.getName().split('-')[0]
            frac_area.setName(domain+"-area_fractions")
            frac_area_type = asearch.child_by_name(frac_area, "field evaluator type")
            if frac_area_type.getValue() == "surface balance area fractions":
                frac_area_type.setValue("area fractions, three components")

    # water source
    try:
        water_source = asearch.find_path(xml, ['state', 'field evaluators', eval_name])
    except aerrors.MissingXMLError:
        pass
    else:
        water_source_type = asearch.child_by_name(water_source, "field evaluator type")
        if water_source_type.getValue() == "surface balance":
            water_source_type.setValue("surface energy balance, three components")

    # subgrid albedos
    for name in ['surface-subgrid_albedos', 'surface_column:*-subgrid_albedos']:
        try:
            albedo = asearch.find_path(xml, ['state', 'field evaluators', name])
        except aerrors.MissingXMLError:
            pass
        else:
            dname = name.split('-')[0]
            albedo.setName(dname+'-albedos')
            albedo_type = asearch.child_by_name(albedo, "field evaluator type")
            if albedo_type.getValue() == "albedo":
                albedo_type.setValue("subgrid albedos, three components")

def create_landcover(xml, eval_name='water_source', water_transition_depth=0.02):
    """Adds a default land-cover section that will be the same as the Arctic defaults."""
    ic_list = asearch.find_path(xml, ['state', 'initial conditions'])
    try:
        lc = asearch.child_by_name(ic_list, "land cover types")
    except aerrors.MissingXMLError:
        pass
    else:
        return # already exists

    # create a single land-cover type on the entire surface domain
    def add_to_lc(eval_name, param_name, param_default, lc, param_name_new=None):
        if param_name_new is None:
            param_name_new = param_name

        # add the dessicated zone thickness
        try:
            evaluator = asearch.find_path(xml, ['state', 'field evaluators', eval_name])
        except aerrors.MissingXMLError:
            pval = None
        else:
            try:
                param = asearch.child_by_name(evaluator, param_name)
            except aerrors.MissingXMLError:
                pval = None
            else:
                pval = param.getValue()
                evaluator.remove(param)

        if pval is None:
            pval = param_default
        lc.append(parameter.DoubleParameter(param_name_new, pval))

    lc = ic_list.sublist('land cover types').sublist('surface domain')
    add_to_lc('water_source', 'dessicated zone thickness [m]', 0.1, lc)
    add_to_lc('water_source', 'roughness length of bare ground [m]', 0.04, lc)
    add_to_lc('water_source', 'roughness length of snow-covered ground [m]', 0.004, lc, 'roughness length of snow [m]')
    add_to_lc('water_source', 'snow-ground transitional depth [m]', 0.02, lc, 'snow transition depth [m]')
    add_to_lc('water_source', 'water-ground transitional depth [m]', water_transition_depth, lc, 'water transition depth [m]')
    add_to_lc('surface-subgrid_albedos', 'albedo ground surface [-]', 0.135, lc, 'albedo of bare ground [-]')
    add_to_lc('surface-subgrid_albedos', 'emissivity tundra [-]', 0.92, lc, 'emissivity of bare ground [-]')


def arctic_seb(xml, seb_new=False):
    if seb_new:
        for name in ['snow-source_sink', 'snow_column:*-source_sink', 'water_source', 'column:*-water_source']:
            try:
                # check for water source
                asearch.find_path(xml, ['state', 'field evaluators', name])
            except aerrors.MissingXMLError:
                pass
            else:
                seb_threecomponent(xml, name)
                create_landcover(xml, name, 0.002)
    else:
        for name in ['snow-source_sink', 'snow_column:*-source_sink', 'water_source', 'column:*-water_source']:
            try:
                # check for water source
                asearch.find_path(xml, ['state', 'field evaluators', name])
            except aerrors.MissingXMLError:
                pass
            else:
                seb_twocomponent(xml, name)
                create_landcover(xml, name, 0.02)


    
def priestley_taylor(xml):
    eval_list = asearch.find_path(xml, ['state', 'field evaluators'])
    for ev in eval_list:
        ev_type = asearch.child_by_name(ev, 'field evaluator type')
        if ev_type.getValue() == 'potential evapotranspiration':
            ev_type.setValue('potential evapotranspiration, Priestley-Taylor')

            # net radiation is the new key, old used shortwave by default
            default = True
            try:
                sw_key = asearch.child_by_name(ev, 'shortwave radiation key')
            except aerrors.MissingXMLError:
                pass
            else:
                sw_key.setName('net radiation key')
                default = False

            try:
                sw_suffix = asearch.child_by_name(ev, 'shortwave radiation key suffix')
            except aerrors.MissingXMLError:
                pass
            else:
                sw_key.setName('net radiation key suffix')
                default = False

            if default:
                ev.append(parameter.StringParameter('net radiation key suffix', 'shortwave_radiation'))

        if ev.getName() == 'surface-air_temperature_inter':
            # this name is dead, now just surface-temperature (which is what it was intended to be physically)
            ev.setName('surface-temperature')
            

def update(xml, seb_new=False):
    alias_mesh_target(xml)
    arctic_seb(xml, seb_new)
    priestley_taylor(xml)
    

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Fix a number of changes from ATS input spec 1.1 to release 1.2")
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
    

    
