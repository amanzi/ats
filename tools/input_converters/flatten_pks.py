"""Flattens PKs, creating cycle driver list.

Part of 0.86 to 0.9X conversion.
"""

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
from amanzi_xml.common.parameter_list import ParameterList
import amanzi_xml.common.parameter as parameter
from amanzi_xml.utils import errors as aerrors


def coordinator_to_cycle_driver(xml):
    cycle_driver = asearch.child_by_name(xml, "coordinator")
    cycle_driver.set("name", "cycle driver")
    cycle_driver_pks = ParameterList("PK tree")
    cycle_driver.append(cycle_driver_pks)

    return cycle_driver_pks

def flatten(pks, flat_pks, cd_pks):
    while len(pks) > 0:
        pk_name = pks.getchildren()[0].get("name")
        new_cd = ParameterList(pk_name)
        pk = pks.pop(pk_name)
        flat_pks.append(pk)
        print(f'Flattened PK: {pk_name}')
        new_cd.append(parameter.StringParameter("PK type", asearch.child_by_name(pk, "PK type").get("value")))
        cd_pks.append(new_cd)

        try:
            subpks = pk.pop("PKs")
        except aerrors.MissingXMLError:
            pass
        else:
            flatten(subpks, flat_pks, new_cd)

def flatten_pks(xml):
    try:
        cycle_driver = asearch.child_by_name(xml, "cycle driver")
    except aerrors.MissingXMLError:
        pass
    else:
        return
    
    pks = asearch.child_by_name(xml, "PKs")
    cd_pks = coordinator_to_cycle_driver(xml)
    flat_pks = ParameterList("PKs")
    flatten(pks, flat_pks, cd_pks)
    pks.extend(flat_pks)
        


if __name__ == "__main__":
    if "-h" in sys.argv or "--help" in sys.argv or "--h" in sys.argv:
        print("Usage: python flatten_pks.py INFILE OUTFILE")
        sys.exit(0)

    outfile = sys.argv.pop(-1)
    infile = sys.argv.pop(-1)
    if (infile == "-i"):
        infile = outfile

    xml = aio.fromFile(infile)
    flatten_pks(xml)
    aio.toFile(xml, outfile)
    sys.exit(0)
        
    
    
    





