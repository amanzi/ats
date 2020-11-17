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


def checkManning(xml):
    fe_list = asearch.gen_by_path(xml, ["state","field evaluators"])
    for eval in fe_list:
        if eval.get("name").endswith('manning_coefficient'):
            eval_type = eval.getElement('field evaluator type')
            if eval_type.getValue() == 'independent variable':
                func_reg = eval.getElement("function")
                for reg in func_reg:
                    comp_entries = asearch.findall_name(reg, ['components'])
                    if len(comp_entries) > 1:
                        # previous iterations of this script were broken...
                        for entry in comp_entries[1:]:
                            reg.remove(entry)
                        

                    fixed = False
                    if not fixed:
                        try:
                            comp = reg.getElement('component')
                        except aerrors.MissingXMLError:
                            pass
                        else:
                            reg.pop('component')
                            assert not any(el.getName() == 'component' for el in reg)
                            reg.append(parameter.ArrayStringParameter('components', ['cell', 'boundary_face']))
                            fixed = True

                    if not fixed:
                        try:
                            comp = reg.getElement('components')
                        except aerrors.MissingXMLError:
                            pass
                        else:
                            reg.pop('components')
                            assert not any(el.getName() == 'components' for el in reg)
                            reg.append(parameter.ArrayStringParameter('components', ['cell', 'boundary_face']))
                            fixed = True

                    if not fixed:
                        print(reg.__str__().decode('utf-8'))
                        raise aerrors.MissingXMLError('Missing "component" or "components"')


def update(xml):
    checkManning(xml)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Fix a number of changes from ATS input spec 0.88 to 1.0")
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
    

    
