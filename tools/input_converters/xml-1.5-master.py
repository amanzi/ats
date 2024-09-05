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


def enforceDtHistory(xml):
    """Find and revert the timestep from file option, moving it to the cycle driver list."""
    ti_file_pars = None
    for ti in asearch.findall_name(xml, "time integrator"):
        if asearch.child_by_name(ti, "timestep controller type").getValue() == "from file":
            for child in ti:
                if child.getName().startswith("timestep controller") and child.getName() != "timestep controller from file parameters":
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


def update(xml):
    enforceDtHistory(xml)


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
