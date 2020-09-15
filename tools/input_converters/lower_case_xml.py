"""
Converts native-spec xml files from old capitalized form to new lower-case form.

Usage:

export AMANZI_SRC_DIR=/path/to/my/amanzi/repo
python lower_case_xml.py my_native_file.xml ... my_other_native_file.xml

NOTE: all work is done in-place!
"""


import sys,os
sys.path.append(os.path.join(os.environ["AMANZI_SRC_DIR"],'tools', 'amanzi_xml'))

import amanzi_xml.utils.errors as aerrors
import amanzi_xml.utils.search as asearch
import amanzi_xml.utils.io as xml_io




changes_regions = ["Coordinate",
                   "Low Coordinate",
                   "High Coordinate",
                   "Location",
                   "Direction",
                   "Tolerance",
                   "Entity",
                   "File",
                   "Format",
                   "Label",
                   "Value",
                   "Number of points",
                   "Points",
                   "Union",
                   "Entity GIDs",
                   "Region: Box",
                   "Region: Plane",
                   "Region: Polygon",
                   "Region: Labeled Set",
                   "Region: Color Function",
                   "Region: Point",
                   "Region: Logical",
                   "Region: Enumerated Set",
                   "Region: All",
                   "Region: Boundary",
                   "Region: Box Volume Fractions"
                   ]

changes_mesh = ["Number of Cells",
                "Domain Low Corner",
                "Domain High Corner",
                "Domain Low Coordinate",
                "Domain High Coordinate",
                "Framework",
                "Verify Mesh",
                "File",
                "Format",
                "Expert",
                "Read Mesh File",
                "Generate Mesh",
                "Surface Mesh"
                ]


def lower_case(xml):
    """Converts an xml object, in-place"""

    try:
        regions = asearch.child_by_name(xml, "Regions")
    except aerrors.MissingXMLError:
        pass

    else:
        for label in changes_regions:
            asearch.change_name(regions, label, label.lower(), allow_multiple=True)

        asearch.change_name(xml, ["Region: Plane","location"], "point", allow_multiple=True)
        asearch.change_name(xml, ["Region: Plane","direction"], "normal", allow_multiple=True)
        regions.setName("regions")
        
    try:
        mesh = asearch.child_by_name(xml, "Mesh")
    except aerrors.MissingXMLError:
        pass
    else:
        for label in changes_mesh:
            asearch.change_name(mesh, label, label.lower(), True)
        mesh.setName("mesh")

    try:
        domain = asearch.child_by_name(xml, "Domain")
    except aerrors.MissingXMLError:
        pass
    else:
        xml.remove(domain)


if __name__ == "__main__":
    fnames = [f for f in sys.argv if f.endswith('.xml')]
    for f in fnames:
        try:
            xml = xml_io.fromFile(f)
        except aerrors.NotNativeSpecError:
            pass
        except RuntimeError,msg:
            print "Attempted to convert invalid xml file:", f
            print "  Error:", msg

        else:        
            lower_case(xml)
            xml_io.toFile(xml, f)





