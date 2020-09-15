#!/usr/bin/env python

"""ATS input converter from 0.87 to dev"""

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
from amanzi_xml.common import parameter, parameter_list

def get_gmres(error_tolerance=1.e-6):
    plist = parameter_list.ParameterList("gmres parameters")
    plist.setParameter("preconditioning strategy", "string", "left")
    plist.setParameter("error tolerance", "double", error_tolerance)
    plist.setParameter("convergence criteria", "Array(string)", ["relative residual", "make one iteration"])
    plist.setParameter("maximum number of iteration", "int", 80)
    return plist

def update(xml, clobber=True):
    for lin_op_list in asearch.findall_name(xml, "linear solver"):
        if lin_op_list.getElement("iterative method").getValue() == "nka":
            lin_op_list.setParameter("iterative method", "string", "gmres")
            lin_op_list.pop("nka parameters")
            lin_op_list.append(get_gmres())
        elif clobber and lin_op_list.getElement("iterative method").getValue() == "gmres":
            lin_op_list.pop("gmres parameters")
            lin_op_list.append(get_gmres())
            


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Update NKA-based linear operator lists to GMRES based lists.  Note this is not valid for old specs (requires >= 0.87)")
    parser.add_argument("infile", help="input filename")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-i", "--inplace", action="store_true", help="fix file in place")
    group.add_argument("-o", "--outfile", help="output filename")

    args = parser.parse_args()

    print "Converting file: %s"%args.infile
    xml = aio.fromFile(args.infile, True)
    update(xml)
    if args.inplace:
        aio.toFile(xml, args.infile)
    else:
        aio.toFile(xml, args.outfile)
    sys.exit(0)
    

    
