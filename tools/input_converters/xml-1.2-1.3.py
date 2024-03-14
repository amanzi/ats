#!/usr/bin/env python3
"""ATS input converter from version 1.2 to version 1.3"""

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

def manning_boundary_face(xml):
    """Converts manning coeficients from boundary face and cell back to just cell"""
    matches = asearch.findall_path(xml, ["state", "field evaluators", "surface-manning_coefficient", "components"]) \
        + asearch.findall_path(xml, ["state", "field evaluators", "surface_star-manning_coefficient", "components"])
    for mann in matches:
        mann.setValue(["cell"])

def mobile_depth(xml):
    """Overland conductivity evals now use 'mobile depth key' instead of 'depth key'"""
    matches = asearch.findall_path(xml, ["state", "field evaluators", "surface-overland_conductivity", "depth key"]) \
        + asearch.findall_path(xml, ["state", "field evaluators", "surface_star-overland_conductivity", "depth key"])
    for md in matches:
        md.setName("mobile depth key")

def mass_flux_water_flux(xml):
    """Mass flux --> water flux to be consistent with sources, distinguish for reactive transport (what mass, C or H2O?)"""
    # boundary conditions in flow pks
    for match in asearch.findall_path(xml, ["boundary conditions", "mass flux",]):
        match.setName("water flux")
    for match in asearch.findall_path(xml, ["boundary conditions", "outward mass flux",]):
        match.setName("outward water flux")

    # transport flux
    try:
        field_evals = asearch.find_path(xml, ["state", "field evaluators"])
    except aerrors.MissingXMLError:
        pass
    else:
        for e in field_evals:
            if "mass_flux" in e.getName():
                e.setName(e.getName().replace("mass_flux","water_flux"))

    for match in asearch.findall_path(xml, ["PKs", "darcy flux key"]):
        # richards PK key for flux
        if "mass_flux" in match.getValue():
            match.setValue(match.getValue().replace("mass_flux", "water_flux"))

    # mass flux key --> water flux key in transport and energy
    for match in asearch.findall_path(xml, ["PKs", "mass flux key"]):
        match.setName("water flux key")
        if match.getValue() == "mass_flux":
            match.setValue("water_flux")
        elif match.getValue() == "surface-mass_flux":
            match.setValue("surface-water_flux")
    for match in asearch.findall_path(xml, ["PKs", "mass flux key suffix"]):
        match.setName("water flux key suffix")
        if match.getValue() == "mass_flux":
            match.setValue("water_flux")
        elif match.getValue() == "surface-mass_flux":
            match.setValue("surface-water_flux")

    # observations of mass flux are common
    for match in asearch.findall_path(xml, ["observations", "variable"]):
        if "mass_flux" in match.getValue():
            match.setValue(match.getValue().replace("mass_flux", "water_flux"))
            
def update(xml):
    """generic update calls all needed things""" 
    manning_boundary_face(xml)
    mobile_depth(xml)
    mass_flux_water_flux(xml)
    

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
    

    
