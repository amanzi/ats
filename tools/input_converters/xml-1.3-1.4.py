#!/usr/bin/env python3
"""ATS input converter from version 1.3 to version 1.4"""

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

def new_state(xml):
    """Converts manning coeficients from boundary face and cell back to just cell"""
    try:
        state_evals = asearch.find_path(xml, ["state", "field evaluators"], no_skip=True)
    except aerrors.MissingXMLError:
        try:
            state_evals = asearch.find_path(xml, ["state", "evaluators"], no_skip=True)
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

    for eval_deps in asearch.findall_path(state_evals, ["evaluator dependencies"]):
        eval_deps.setName("dependencies")

def pk_initial_timestep(xml, has_orig=False):
    """'initial time step' --> 'initial time step [s]'"""

    def get_this_dt(pk_tree, pks):
        """Finds dt in THIS list"""
        pk_list = asearch.child_by_name(pks, pk_tree.getName())
        init_ts = list(asearch.children_by_name(pk_list, "initial time step"))
        init_ts_s = list(asearch.children_by_name(pk_list, "initial time step [s]"))
        dt = -1
        if len(init_ts_s) == 1:
            dt = init_ts_s[0].getValue()
            pk_list.remove(init_ts_s[0])
            if len(init_ts) == 1:
                pk_list.remove(init_ts[0])
        elif len(init_ts) == 1:
            dt = init_ts[0].getValue()
            pk_list.remove(init_ts[0])
        return dt

    def get_dt(pk_tree, pks):
        """Finds dt in any PK, if it exists, through recursion of children or local plist"""
        dt = get_this_dt(pk_tree, pks)
        if len(pk_tree) == 1:
            # it is a leaf
            return dt
        elif dt > 0:
            # not a leaf, but has a dt
            return dt
        else:
            # not a leaf, recurse
            dts = [get_dt(child, pks) for child in pk_tree if isinstance(child, parameter_list.ParameterList)]
            if any(dt > 0 for dt in dts):
                return min(dt for dt in dts if dt > 0)
            else:
                return -1

    def generate_with_timestepper(pk_tree, pks):
        """Generator over all PKs that have a time integrator list"""
        for pk in pk_tree:
            if isinstance(pk, parameter_list.ParameterList):
                pk_list = asearch.child_by_name(pks, pk.getName())
                try:
                    ti = asearch.child_by_name(pk_list, "time integrator")
                except aerrors.MissingXMLError:
                    for pk_with_ti in generate_with_timestepper(pk, pks):
                        yield pk_with_ti
                else:
                    yield pk

    pk_tree = asearch.find_path(xml, ["cycle driver", "PK tree"], no_skip=True)
    pks = asearch.find_path(xml, ['PKs',])

    for pk in generate_with_timestepper(pk_tree, pks):
        dt = get_dt(pk, pks)
        if dt > 0:
            pk_list = asearch.child_by_name(pks, pk.getName())
            ti = asearch.child_by_name(pk_list, "time integrator")
            ts_cont_type = asearch.child_by_name(ti, "timestep controller type").getValue()

            ts_cont_pars = None
            if has_orig and ts_cont_type == "from file":
                # this has an orig file, and this is from file, so
                # really we need to put the init dt in the old,
                # original list.
                for child in ti:
                    child_name = child.getName()
                    if child_name.startswith("timestep controller") \
                       and child_name.endswith("parameters") \
                       and child_name != "timestep controller from file parameters":
                        ts_cont_pars = child
            if ts_cont_pars is None:
                ts_cont_pars = asearch.child_by_name(ti, "timestep controller "+ts_cont_type.strip()+" parameters")

            try:
                dt = asearch.child_by_name(ts_cont_pars, "initial time step [s]")
            except aerrors.MissingXMLError:
                ts_cont_pars.setParameter("initial time step [s]", 'double', float(dt))

def top_cell_evals(xml):
    evals = asearch.find_path(xml, ['state', 'evaluators'], no_skip=True)
    for ev in evals:
        ev_type = asearch.child_by_name(ev, 'evaluator type')
        if ev_type.getValue() == 'surface top cell evaluator':
            ev_type.setValue('surface from top cell evaluator')


def density_units(xml):
    """[kg/m^3] --> [kg m^-3] to be consistent with all other units in ATS"""
    for d in asearch.findall_path(xml, ["molar mass [kg/mol]"]):
        d.setName("molar mass [kg mol^-1]")
    for d in asearch.findall_path(xml, ["molar mass [g/mol]"]):
        d.setName("molar mass [g mol^-1]")
    for d in asearch.findall_path(xml, ["density [kg/m^3]"]):
        d.setName("density [kg m^-3]")
    for d in asearch.findall_path(xml, ["density [mol/m^3]"]):
        d.setName("density [mol m^-3]")


def rh_to_vp(xml):
    """Converts relative humidity to vapor pressure."""
    import logging

    for ev in asearch.findall_path(xml, ["state", "evaluators", "surface-relative_humidity"], no_skip=True):
        ev.setName("surface-vapor_pressure_air")
        for fname, hname in zip(["function-tabular", "function-bilinear-and-time"], ["y header", "value header"]):
            try:
                y_header = asearch.find_path(ev, ["function", "domain", "function", fname, hname], no_skip=True)
                y_header.setValue("vapor pressure air [Pa]")
            except aerrors.MissingXMLError:
                logging.warning("Manually update is required! Make sure to use the name of the air vapor pressure in your forcing dataset for the \"value header\" of \"surface-vapor_pressure_air\".")
            else:
                logging.warning("Changing relative_humidity --> vapor_pressure; please update your forcing data to include vapor pressure rather than relative humidity.  One way to do this is to run `$ATS_SRC_DIR/tools/utils/rh_to_vp.py --inplace path/to/daymet.h5`")
    
    for match in asearch.findall_path(xml, ["observations", "variable"]):
        if "relative_humidity" in match.getValue():
            match.setValue(match.getValue().replace("relative_humidity", "vapor_pressure_air"))
            logging.warning("The surface-relative_humidity in observations has been changed to surface-vapor_pressure_air. You may want to change the corresponding output variable name, e.g., from \"relative humidity [-]\" to \"vapor pressure air [Pa]\".")

    # delete relative humidity key if exists
    try:
        _ = asearch.remove_element(xml, "relative humidity key", allow_multiple=True)
    except aerrors.MissingXMLError:
        pass


def pk_flow_reactive_transport(xml):
    pk_tree = asearch.find_path(xml, ["cycle driver", "PK tree"], no_skip=True)
    for pk_type_in_tree in asearch.findall_path(pk_tree, ["PK type"]):
        if pk_type_in_tree.getValue() == "flow reactive transport":
            pk_type_in_tree.setValue("subcycling MPC")

    for pk in asearch.find_path(xml, ["PKs"], no_skip=True):
        pk_type = asearch.child_by_name(pk, "PK type")
        if pk_type.getValue() == "flow reactive transport":
            pk_type.setValue("subcycling MPC")
            try:
                pk.pop("master PK index")
            except asearch.MissingXMLError:
                pass
            else:
                pk.setParameter("subcycle", "Array(int)", [0,1])

    for pk in asearch.find_path(xml, ["PKs"], no_skip=True):
        pk_type = asearch.child_by_name(pk, "PK type")
        if pk_type.getValue() == "transport ATS":
            try:
                tests = asearch.child_by_name(pk, "enable internal tests")
            except aerrors.MissingXMLError:
                pass
            else:
                if tests.getValue() == "yes":
                    pk.pop("enable internal tests")
                    pk.setParameter("enable internal tests", "bool", True)
                elif tests.getValue() == "no":
                    pk.pop("enable internal tests")
                    pk.setParameter("enable internal tests", "bool", False)
            
                
def update(xml, has_orig=False):
    """generic update calls all needed things""" 
    new_state(xml)
    pk_initial_timestep(xml, has_orig)
    top_cell_evals(xml)
    density_units(xml)
    rh_to_vp(xml)
    pk_flow_reactive_transport(xml)
    
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("infile", help="input filename")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-i", "--inplace", action="store_true", help="fix file in place")
    group.add_argument("-o", "--outfile", help="output filename")

    args = parser.parse_args()

    # check for orig file
    has_orig = False
    if args.infile.endswith('.xml'):
        orig_filename = args.infile[:-4]+"_orig.xml"
        if os.path.isfile(orig_filename):
            has_orig = True

    print("Converting file: %s"%args.infile)
    xml = aio.fromFile(args.infile, True)
    update(xml)
    if args.inplace:
        aio.toFile(xml, args.infile)
    else:
        aio.toFile(xml, args.outfile)
    sys.exit(0)
    

    
