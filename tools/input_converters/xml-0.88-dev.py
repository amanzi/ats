#!/usr/bin/env python3
"""ATS input converter from 0.88 to master"""

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

def change_name(xml, old, new):
    try:
        asearch.change_name(xml, old, new)
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

def checkManning(xml):
    fe_list = asearch.find_path(xml, ["state","field evaluators"])
    for eval in fe_list:
        if eval.getName().endswith('manning_coefficient'):
            eval_type = eval.getElement('field evaluator type')
            if eval_type.getValue() == 'independent variable':
                func_reg = eval.getElement("function")
                for reg in func_reg:
                    fixed = False
                    if not fixed:
                        try:
                            comp = reg.getElement('component')
                        except aerrors.MissingXMLError:
                            pass
                        else:
                            reg.pop('component')
                            reg.append(parameter.ArrayStringParameter('components', ['cell', 'boundary_face']))
                            fixed = True

                    if not fixed:
                        try:
                            comp = reg.getElement('components')
                        except aerrors.MissingXMLError:
                            pass
                        else:
                            reg.pop('components')
                            reg.append(parameter.ArrayStringParameter('components', ['cell', 'boundary_face']))
                            fixed = True

                    if not fixed:
                        print(reg.__str__().decode('utf-8'))
                        raise aerrors.MissingXMLError('Missing "component" or "components"')


def fixSnow(xml):
    # todo: remove "include dt factor"
    # todo: rename "dt factor" --> "dt factor [s]"
    # todo: add "swe density factor [-]" (default is 10)
    # todo: rename "include density factor" --> "include density"
    raise NotImplementedError("fix snow")

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
           

def update(xml):
    vanGenuchtenParams(xml)
    checkManning(xml)
    # NOTE: these will get added when subgrid pull request is done
    #fixSnow(xml)
    #fixSubgrid(xml)
    mergePreconditionerLinearSolver(xml)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Fix a number of changes from ATS input spec 0.88 to master")
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
    

    
