#!/usr/bin/env python3
"""ATS input converter from main development branch to make tests that are friendlier for testing tpetra branch."""

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


def make_constant(fe_list, var, val, alternate=None):
    # already done
    if fe_list.isElement(var) and fe_list.sublist(var).getElement("evaluator type").getValue() == "independent variable constant":
        return
    # not there, no alternate name
    if not fe_list.isElement(var) and alternate is None:
        return
    # not there, alternate not there
    if not fe_list.isElement(var) and alternate is not None and not fe_list.isElement(alternate):
        return
    
    if fe_list.isElement(var):
        fe_list.pop(var)
    elist = fe_list.sublist(var)
    elist.setParameter("evaluator type", "string", "independent variable constant")
    elist.setParameter("value", "double", val)


def mfd(xml):
    for par in asearch.findall_path(xml, ["PKs", "discretization primary"]):
        if par.getValue() == "fv: default" or par.getValue() == "mfd: two-point flux approximation":
            pass
        else:
            par.setValue("mfd: default")


def update(xml):
    """generic update calls all needed things"""
    fe_list = asearch.find_path(xml, ["state", "evaluators"], no_skip=True)
    make_constant(fe_list, 'viscosity_liquid', 8.9e-4)
    make_constant(fe_list, 'molar_density_liquid', 55500, 'mass_density_liquid')
    make_constant(fe_list, 'mass_density_liquid', 1000, 'molar_density_liquid')
    make_constant(fe_list, 'surface-mass_density_liquid', 900, 'surface-molar_density_liquid')
    make_constant(fe_list, 'surface-molar_density_liquid', 55500, 'surface-mass_density_liquid')
    make_constant(fe_list, 'surface-source_molar_density', 55500)
    mfd(xml)


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
    

    
