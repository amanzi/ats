#!/usr/bin/env python

import sys,os
sys.path.append(os.path.join(os.environ['AMANZI_SRC_DIR'], 'tools', 'amanzi_xml'))
from amanzi_xml.utils import errors
from amanzi_xml.utils import io
from amanzi_xml.utils import search as asearch
from amanzi_xml.common import parameter_list, parameter


def removeID(xml):
    try:
        xml.attrib.pop('id')
    except KeyError:
        pass
    for el in xml:
        if 'type' in el.keys() and el.attrib['type'] == 'ParameterList':
            removeID(el)
        else:
            try:
                el.attrib.pop('id')
            except KeyError:
                pass


def fixPorosity(xml):
    """Fix the compressible porosity list"""
    try:
        asearch.find_by_path(xml, ["Main","Regions","computation domain moss"])
    except errors.MissingXMLError:
        moss = False
    else:
        moss = True

    try:
        poro = asearch.find_by_path(xml, ["Main","state","field evaluators","porosity"])
    except errors.MissingXMLError:
        print "no porosity evaluator"
        return

    try:
        params = asearch.find_by_path(poro, ["compressible porosity model parameters",])
    except errors.MissingXMLError:
        print "porosity is not compressible"
        return

    if params.getchildren()[0].get("type") == "ParameterList":
        print "porosity has already been fixed"
        return

    comp_el = params.getchildren().pop()
    assert len(params.getchildren()) == 0

    if moss:
        moss_list = parameter_list.ParameterList("moss")
        moss_list.setParameter("region", "string", "computational domain moss")
        moss_list.setParameter("pore compressibility", "double", XXX)
        params.getchildren().append(moss_list)

        peat_list = parameter_list.ParameterList("peat")
        peat_list.setParameter("region", "string", "computational domain peat")
        peat_list.setParameter("pore compressibility", "double", XXX)
        params.getchildren().append(peat_list)

        up_min_list = parameter_list.ParameterList("upper mineral")
        up_min_list.setParameter("region", "string", "computational domain upper mineral")
        up_min_list.setParameter("pore compressibility", "double", XXX)
        params.getchildren().append(up_min_list)

        low_org_list = parameter_list.ParameterList("lower organic")
        low_org_list.setParameter("region", "string", "computational domain lower organic")
        low_org_list.setParameter("pore compressibility", "double", XXX)
        params.getchildren().append(low_org_list)

        low_min_list = parameter_list.ParameterList("lower mineral")
        low_min_list.setParameter("region", "string", "computational domain lower mineral")
        low_min_list.setParameter("pore compressibility", "double", XXX)
        params.getchildren().append(low_min_list)

        deep_min_list = parameter_list.ParameterList("deep mineral")
        deep_min_list.setParameter("region", "string", "computational domain deep mineral")
        deep_min_list.setParameter("pore compressibility", "double", XXX)
        params.getchildren().append(deep_min_list)

    else:
        cp_list = parameter_list.ParameterList("computational domain")
        cp_list.setParameter("region", "string", "computational domain")
        cp_list.setParameter("pore compressibility", "double", comp_el.get("value"))
        params.getchildren().append(cp_list)


def migrate(new, old, name, oldname=None):
    """move a param from old to new, if it exists"""
    if oldname is None:
        oldname = name
    try:
        el = old.pop(oldname)
    except errors.MissingXMLError:
        pass
    else:
        el.set("name", name)
        new.getchildren().append(el)


def _fixTIList(ti_list):
    if ti_list.isElement("solver type"):
        # already new style
        return

    old_pars = parameter_list.ParameterList("time integrator")
    while len(ti_list.getchildren()) > 0:
        old_pars.getchildren().append(ti_list.getchildren().pop())

    # ti
    migrate(ti_list, old_pars, "max preconditioner lag iterations")
    migrate(ti_list, old_pars, "extrapolate initial guess")

    # nonlinear solver
    nl_solver = old_pars.pop("nonlinear solver").get("value")

    if nl_solver == "NKA":
        ti_list.setParameter("solver type", "string", "nka")
        nl_plist = ti_list.sublist("nka parameters")
        migrate(nl_plist, old_pars, "max du growth factor")
        migrate(nl_plist, old_pars, "max error growth factor")
        migrate(nl_plist, old_pars, "max divergent iterations")
        migrate(nl_plist, old_pars, "lag iterations")
        nl_plist.setParameter("modify correction", "bool", True)
        migrate(nl_plist, old_pars, "monitor")
        migrate(nl_plist, old_pars, "monitor", "convergence monitor")
        nl_plist.setParameter("monitor", "string", "monitor residual")


    elif nl_solver == "NKA BT":
        ti_list.setParameter("solver type", "string", "nka_bt_ats")
        nl_plist = ti_list.sublist("nka_bt_ats parameters")
        migrate(nl_plist, old_pars, "nka lag iterations")
        migrate(nl_plist, old_pars, "max backtrack steps")
        migrate(nl_plist, old_pars, "max total backtrack steps")
        migrate(nl_plist, old_pars, "backtrack lag")
        migrate(nl_plist, old_pars, "last backtrack iteration")
        migrate(nl_plist, old_pars, "backtrack factor")
        migrate(nl_plist, old_pars, "backtrack tolerance")

    migrate(nl_plist, old_pars, "nonlinear tolerance")
    migrate(nl_plist, old_pars, "diverged tolerance")
    migrate(nl_plist, old_pars, "limit iterations")
    nl_plist.sublist("VerboseObject").setParameter("Verbosity Level", "string",
                                                   old_pars.sublist("VerboseObject").getElement("Verbosity Level").get("value"))

    migrate(ti_list, old_pars, "VerboseObject")

    # timestep controller
    ts_ctrl = old_pars.pop("timestep controller type").get("value")
    assert ts_ctrl == "smarter"
    ti_list.setParameter("timestep controller type", "string", ts_ctrl)
    ts_list = ti_list.sublist("timestep controller smarter parameters")
    migrate(ts_list, old_pars, "max iterations")
    migrate(ts_list, old_pars, "min iterations")
    migrate(ts_list, old_pars, "time step reduction factor")
    migrate(ts_list, old_pars, "time step increase factor")
    migrate(ts_list, old_pars, "max time step increase factor")
    migrate(ts_list, old_pars, "max time step")
    migrate(ts_list, old_pars, "min time step")
    migrate(ts_list, old_pars, "growth wait after fail")
    migrate(ts_list, old_pars, "count before increasing increase factor")

def fixTIList(xml):
    for ti_list in asearch.findall_name(xml, ["time integrator",]):
        _fixTIList(ti_list)


def _fixPC(pc, pctype, force=False):
    if (not pc.isElement("preconditioner type")) or force:
        pc.setParameter("preconditioner type","string",pctype)

def _fixLinOp(itr):
    if not itr.isElement("iterative method"):
        itr.setParameter("iterative method", "string", "nka")
    else:
        itm = itr.getElement("iterative method").get("value")
        if itm == "nka":
            lo = itr.sublist("nka parameters")
        elif itm == "gmres":
            lo = itr.sublist("gmres parameters")
        else:
            raise RuntimeError("unknown lin op %s"%itm)

        if not lo.isElement("error tolerance"):
            lo.setParameter("error tolerance", "double", 1.e-6)
        if not lo.isElement("maximum number of iterations"):
            lo.setParameter("maximum number of iterations", "int", 20)

def _fixDiffusion(diff, req_pc, req_face):
    if req_face:
        face = diff.sublist("consistent face solver")
        _fixLinOp(face)
        _fixPC(face.sublist("preconditioner"), "block ilu")

    if req_pc:
        _fixPC(diff.sublist("preconditioner"), "boomer amg")
    else:
        try: diff.pop("preconditioner")
        except errors.MissingXMLError: pass

    if diff.getElement("MFD method").get("value") == "two point flux":
        diff.setParameter("MFD method", "string", "two point flux approximation")

def fixSolvers(xml):
    for diff in asearch.findall_name(xml, "Diffusion"):
        _fixDiffusion(diff, False, True)

    need_pc = True
    for pk in asearch.findall_name(xml, "PK type"):
        if pk.get("value") == "new permafrost model no SC" or \
           pk.get("value") == "subsurface permafrost" or \
           pk.get("value") == "coupled water":
            need_pc = False
    for diff in asearch.findall_name(xml, "Diffusion PC"):
        _fixDiffusion(diff, need_pc, True)

    for solver in asearch.findall_name(xml, "Coupled PC"):
        _fixPC(solver.sublist("preconditioner"), "boomer amg", True)

    for solver in asearch.findall_name(xml, "Coupled Solver"):
        _fixLinOp(solver)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Fix a number of changes from ATS input spec 0.83 to 0.86")
    parser.add_argument("infile", help="input filename")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-i", "--inplace", action="store_true", help="fix file in place")
    group.add_argument("-o", "--outfile", help="output filename")

    args = parser.parse_args()

    # read the infile
    xml = io.fromFile(args.infile)

    # fix
    removeID(xml)
    fixPorosity(xml)
    fixTIList(xml)
    fixSolvers(xml)

    # write the outfile
    if args.inplace:
        io.toFile(xml, args.infile)
    else:
        io.toFile(xml, args.outfile)
    sys.exit(0)

