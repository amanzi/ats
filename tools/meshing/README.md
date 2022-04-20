Meshing Workflows for ATS Mesh Development
==============================================

There are currently four frequently used approaches for generating 3D ATS meshes.

- `extrude`
- `meshing_ats`
- GIS/python-based workflow https://github.com/ecoon/watershed-workflow
- LaGrit/TetGen/others

Anything that can write an Exodus II file or can be converted to a
valid unstructured Exodus II file (e.g. VTK) is usable.  These are
just our most commonly used tools.


extrude
----------

`extrude` is a C++ code that can handle large numbers of triangles
quite easily.  It assumes you have a triangulated surface mesh, and
provides tools for extruding (in the vertical) to generate triangular
primatic meshes.

See `${ATS_SRC_DIR}/tools/meshing_ats/extrude`


meshing_ats
-------------

`meshing_ats` is a python-based code that can handle arbitrary
polyhedra, and is the best choice for making simple (column-based 1D
or 2D transects/hillslopes) meshes.  It is also workable for 3D
meshes, and overlaps with `extrude` in its capabilities there.
It is less performant but more general than extrude.

See `${ATS_SRC_DIR}/tools/meshing_ats/meshing_ats/meshing_ats.py`


GIS/python based workflow
----------------------------

This is the best choice for generating surface meshes for passing off
to `extrude` and/or `meshing_ats` for extrusion.  It can handle DEMs
for you, and can download USGS datasets to do so.  It works with
watershed shapefiles, etc.

See [the repo|https://github.com/ecoon/watershed-workflow.git]


LaGrit/TetGen/others
---------------------

Anything that makes 3D meshes works fine -- we have used LaGrit
extensively for this.

See [the LaGrit Site|https://lagrit.lanl.gov/] and [ATS's project page|https://meshing.lanl.gov/proj/#Arctic]





