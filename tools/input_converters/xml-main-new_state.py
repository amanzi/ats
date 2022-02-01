#!/usr/bin/env python3
"""ATS input converter from main development branch to new-state development branch"""

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

def new_state(xml):
    """Converts manning coeficients from boundary face and cell back to just cell"""
    try:
        state_evals = asearch.find_path(xml, ["state", "field evaluators"])
    except aerrors.MissingXMLError:
        try:
            state_evals = asearch.find_path(xml, ["state", "evaluators"])
        except aerrors.MissingXMLError:
            return
    else:
        state_evals.setName("evaluators")

    for fe in state_evals:
        try:
            fe_type = asearch.child_by_name(fe, "field evaluator type")
        except aerrors.MissingXMLError:
            pass
        else:
            fe_type.setName("evaluator type")
            
def update(xml):
    """generic update calls all needed things""" 
    new_state(xml)

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
    

    
