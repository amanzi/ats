#!/usr/bin/env python3
"""ATS input converter from 1.2 to main development branch"""

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
    matches = asearch.findall_path(xml, ["state", "field evaluators", "surface-manning_coefficient", "components"]) \
        + asearch.findall_path(xml, ["state", "field evaluators", "surface_star-manning_coefficient", "components"])
    for mann in matches:
        mann.setValue(["cell"])

def mobile_depth(xml):
    matches = asearch.findall_path(xml, ["state", "field evaluators", "surface-overland_conductivity", "depth key"]) \
        + asearch.findall_path(xml, ["state", "field evaluators", "surface_star-overland_conductivity", "depth key"])
    for md in matches:
        md.setName("mobile depth key")

def update(xml, seb_new=False):
    manning_boundary_face(xml)
    mobile_depth(xml)
    

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Fix a number of changes from ATS input spec 1.2 to the development branch")
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
    

    
