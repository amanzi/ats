"""Helper function for fixing old style VerboseObject lists"""

import sys, os
try:
    amanzi_xml = os.path.join(os.environ["AMANZI_SRC_DIR"], "tools","amanzi_xml")
except KeyError:
    pass
else:
    if amanzi_xml not in sys.path:
        sys.path.append(amanzi_xml)

from amanzi_xml.utils import search as asearch
from amanzi_xml.utils import errors as aerrors

        
def fixVerboseObject(xml):
    for vo in asearch.generateElementByNamePath(xml, "VerboseObject"):
        vo.set("name", "verbose object")
        try:
            vlev = asearch.childByNamePath(vo, "Verbosity Level")
        except aerrors.MissingXMLError:
            pass
        else:
            vlev.set("name", "verbosity level")

