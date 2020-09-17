"""ATS input converter, moves SEB from the monolithic version in 0.87
and earlier to a newer, modularized version in 0.88."""

import sys, os
try:
    amanzi_xml = os.path.join(os.environ["AMANZI_SRC_DIR"], "tools","amanzi_xml")
except KeyError:
    pass
else:
    if amanzi_xml not in sys.path:
        sys.path.append(amanzi_xml)

import copy

from amanzi_xml.utils import search as asearch
from amanzi_xml.utils import io as aio
from amanzi_xml.utils import errors as aerrors
from amanzi_xml.common import parameter, parameter_list


def add_snow_mesh(xml):
    meshes = asearch.child_by_name(xml, "mesh")
    snow = meshes.sublist("snow")
    snow.append(parameter.StringParameter("mesh type", "aliased"))
    aliased = snow.sublist("aliased parameters")
    aliased.append(parameter.StringParameter("alias", "surface"))

    vis = asearch.child_by_name(xml, "visualization")
    surf_vis = asearch.child_by_name(vis, "surface")
    snow_vis = vis.sublist("snow")
    snow_vis.append(parameter.StringParameter("file name base", "visdump_snow"))
    for p in surf_vis:
        if p.get('name') != "file name base":
            snow_vis.append(copy.copy(p))

def eval_snow_swe():
    swe = parameter_list.ParameterList("snow-swe")
    swe.append(parameter.StringParameter("field evaluator type", "multiplicative evaluator"))
    swe.append(parameter.ArrayStringParameter("evaluator dependencies", ["snow-depth", "snow-density", "snow-cell_volume"]))
    swe.append(parameter.DoubleParameter("coefficient", 1.e-3))
    return swe

def eval_snow_frac_areas():
    fa = parameter_list.ParameterList("surface-fractional_areas")
    fa.append(parameter.StringParameter("field evaluator type", "surface balance area fractions"))
    return fa

def eval_snow_source_sink(dc=None):
    ss = parameter_list.ParameterList("snow-source_sink")
    ss.append(parameter.StringParameter("field evaluator type", "surface balance"))
    ss.append(parameter.BoolParameter("save diagnostic data", True))

    if dc is not None:
        ss.append(dc)
        vo = ss.sublist("verbose object")
        vo.append(parameter.StringParameter("verbosity level", "high"))
    return ss

def eval_longwave():
    lw = parameter_list.ParameterList("surface-incoming_longwave_radiation")
    lw.append(parameter.StringParameter("field evaluator type", "incoming longwave radiation"))
    return lw

def eval_albedo():
    al = parameter_list.ParameterList("surface-subgrid_albedos")
    al.append(parameter.StringParameter("field evaluator type", "albedo"))
    return al


def evals(xml, dc):
    eval_list = asearch.find_by_path(xml, ["state","field evaluators"])
    try:
        lw = asearch.child_by_name(eval_list, "surface-incoming_longwave_radiation")
    except aerrors.MissingXMLError:
        eval_list.append(eval_longwave())

    eval_list.append(eval_snow_frac_areas())
    eval_list.append(eval_snow_swe())
    eval_list.append(eval_snow_source_sink(dc))
    eval_list.append(eval_albedo())

    try:
        snow_precip = asearch.child_by_name(eval_list, "surface-precipitation_snow")
    except aerrors.MissingXMLError:
        pass
    else:
        snow_precip.setName("snow-precipitation")


def pks(xml):
    pk_tree = asearch.child_by_namePath(xml, "cycle driver/PK tree")
    for pk_type in asearch.gen_by_path(pk_tree, "PK type"):
        if pk_type.getValue() == 'surface balance implicit':
            pk_type.setValue('surface balance implicit subgrid')
    
    pks = asearch.child_by_name(xml, "PKs")
    for pk in pks:
        pk_type = asearch.child_by_name(pk, "PK type")
        if pk_type.get('value') == 'permafrost flow':
            try:
                source_term = asearch.child_by_name(pk, "source term")
            except aerrors.MissingXMLError:
                pass
            else:
                source_term.setValue(False)

        elif pk_type.get('value') == 'overland flow with ice':
            try:
                source_term = asearch.child_by_name(pk, "source term")
            except aerrors.MissingXMLError:
                pk.append(parameter.BoolParameter("source term", True))
            else:
                source_term.setValue(True)

            try:
                source_is_diff = asearch.child_by_name(pk, "source term is differentiable")
            except aerrors.MissingXMLError:
                pk.append(parameter.BoolParameter("source term is differentiable", False))
            else:
                source_is_diff.setValue(False)

        elif pk_type.get('value') == 'three-phase energy':
            try:
                source_term = asearch.child_by_name(pk, "source term")
            except aerrors.MissingXMLError:
                pass
            else:
                source_term.setValue(False)

        elif pk_type.get('value') == 'surface energy':
            try:
                source_term = asearch.child_by_name(pk, "source term")
            except aerrors.MissingXMLError:
                pk.append(parameter.BoolParameter("source term", True))
            else:
                source_term.setValue(True)

            try:
                source_is_diff = asearch.child_by_name(pk, "source term is differentiable")
            except aerrors.MissingXMLError:
                pk.append(parameter.BoolParameter("source term is differentiable", True))
            else:
                source_is_diff.setValue(True)

            try:
                source_fd = asearch.child_by_name(pk, "source term finite difference")
            except aerrors.MissingXMLError:
                pk.append(parameter.BoolParameter("source term finite difference", True))
            else:
                source_fd.setValue(True)
                
        elif pk_type.get('value') == 'surface balance implicit':
            pk_seb = parameter_list.ParameterList(pk.get('name'))
            pk_seb.append(parameter.StringParameter('PK type', 'surface balance implicit subgrid'))
            pk_seb.append(parameter.StringParameter('layer name', 'snow'))
            pk_seb.append(parameter.StringParameter('domain name', 'snow'))
            pk_seb.append(parameter.StringParameter('primary variable key', 'snow-depth'))
            pk_seb.append(parameter.StringParameter('conserved quantity key', 'snow-swe'))
            pk_seb.append(parameter.StringParameter('source key', 'snow-source_sink'))
            pk_seb.append(parameter.BoolParameter('source term is differentiable', False))
            pk_seb.append(pk.pop('initial condition'))

            try:
                pk_seb.append(pk.pop('verbose object'))
            except aerrors.MissingXMLError:
                pass

            try:
                dc = pk.pop('debug cells')
            except aerrors.MissingXMLError:
                dc = None
            else:
                pk_seb.append(copy.copy(dc))

            pc = pk_seb.sublist('preconditioner')
            pc.append(parameter.StringParameter('preconditioner type', 'identity'))

            ls = pk_seb.sublist('linear solver')
            ls.append(parameter.StringParameter('iterative method', 'nka'))
            nka = ls.sublist('nka parameters')
            nka.append(parameter.DoubleParameter('error tolerance', 1.e-6))
            nka.append(parameter.IntParameter('maximum number of iterations', 10))
            pks.pop(pk.get('name'))
            pks.append(pk_seb)

    return dc
            

def update_seb(xml):
    add_snow_mesh(xml)
    dc = pks(xml)
    evals(xml, dc)
    


    
