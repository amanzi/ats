#!/usr/bin/env python3
"""ATS input converter to add new observations with signed runoff/on"""

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
import copy


def update(xml):
    observations = asearch.child_by_name(xml, "observations")
    for obs in observations:
        quants = asearch.child_by_name(obs, "observed quantities")

        try:
            asearch.child_by_name(quants, "gross runoff [mol d^-1]")
        except aerrors.MissingXMLError:
            pass
        else:
            continue # did we already do this obs?

        # expected structure?
        if len(quants) < 18:
            # not a water mass balance?
            continue

        # runoff -- add gross runoff and runon
        try:
            runoff = asearch.child_by_name(quants, "net runoff [mol d^-1]")
        except aerrors.MissingXMLError:
            pass
        else:
            runoff_i = next(i for (i,r) in enumerate(quants) if r is runoff)
            
            # gross runoff -- positive outward
            grunoff = copy.deepcopy(runoff)
            grunoff.set("name","gross runoff [mol d^-1]")
            grunoff_mod_list = parameter_list.ParameterList("modifier")
            grunoff_mod_list.setParameter("function type", "string", "standard math")
            grunoff_mod_list.setParameter("operator", "string", "positive")
            grunoff.append(grunoff_mod_list)

            # gross runon -- positive outward
            grunon = copy.deepcopy(runoff)
            grunon.set("name","runon (negative) [mol d^-1]")
            grunon_mod_list = parameter_list.ParameterList("modifier")
            grunon_mod_list.setParameter("function type", "string", "standard math")
            grunon_mod_list.setParameter("operator", "string", "negative")
            grunon.append(grunon_mod_list)

            # insert in reverse order, after the net runoff
            quants.insert(runoff_i+1, grunon)
            quants.insert(runoff_i+1, grunoff)

        # river corridor runoff -- add gross runoff and runon
        try:
            river = asearch.child_by_name(quants, "river discharge [mol d^-1]")
        except aerrors.MissingXMLError:
            pass
        else:
            river_i = next(i for (i,r) in enumerate(quants) if r is river)

            # river-corridor gross runoff -- positive outward
            rgrunoff = copy.deepcopy(river)
            rgrunoff.set("name","river corridor gross runoff [mol d^-1]")
            rgrunoff_mod_list = parameter_list.ParameterList("modifier")
            rgrunoff_mod_list.setParameter("function type", "string", "standard math")
            rgrunoff_mod_list.setParameter("operator", "string", "positive")
            rgrunoff.append(rgrunoff_mod_list)

            # gross runon -- positive outward
            rgrunon = copy.deepcopy(river)
            rgrunon.set("name","river corridor runon (negative) [mol d^-1]")
            rgrunon_mod_list = parameter_list.ParameterList("modifier")
            rgrunon_mod_list.setParameter("function type", "string", "standard math")
            rgrunon_mod_list.setParameter("operator", "string", "negative")
            rgrunon.append(rgrunon_mod_list)

            # insert in reverse order, after the net runoff
            quants.insert(river_i+1, rgrunon)
            quants.insert(river_i+1, rgrunoff)

        # groundwater -- add gross groundwater in and out
        try:
            ground = asearch.child_by_name(quants, "net groundwater flux [mol d^-1]")
        except aerrors.MissingXMLError:
            pass
        else:
            ground_i = next(i for (i,r) in enumerate(quants) if r is ground)

            # gross outward groundwater flux -- positive outward
            ggroundoff = copy.deepcopy(ground)
            ggroundoff.set("name","gross outward groundwater flux [mol d^-1]")
            ggroundoff_mod_list = parameter_list.ParameterList("modifier")
            ggroundoff_mod_list.setParameter("function type", "string", "standard math")
            ggroundoff_mod_list.setParameter("operator", "string", "positive")
            ggroundoff.append(ggroundoff_mod_list)

            # gross inward groundwater flux -- positive outward
            ggroundon = copy.deepcopy(ground)
            ggroundon.set("name","gross inward groundwater flux (negative) [mol d^-1]")
            ggroundon_mod_list = parameter_list.ParameterList("modifier")
            ggroundon_mod_list.setParameter("function type", "string", "standard math")
            ggroundon_mod_list.setParameter("operator", "string", "negative")
            ggroundon.append(ggroundon_mod_list)

            # insert in reverse order, after the net runoff
            quants.insert(ground_i+1, ggroundon)
            quants.insert(ground_i+1, ggroundoff)

        

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
    

    
