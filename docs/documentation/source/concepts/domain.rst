Model Domain
============

Simulations in Amanzi-ATS consist of a coupled set of partial
differential equations, typically solving for the conservation of some
quantity (e.g. water mass or energy), solved on a **domain**.
Mathematically we define a domain to be a submanifold of :math:`R^n`,
where :math:`n \in \{1,2,3\}`.

We use the word domain intentionally to be distinct from the mesh,
which is a discrete representation of the domain.  In Amanzi-ATS,
meshes are always defined in a Cartesian coordinate system, where each
of :math:`\{x,y,z\}` are in units of meters. The domain conceptually
exists without the mesh, though once in the code, the domain is
synonymous with the mesh.  Domains are also used in naming, and
multiple domains may be represented by the same mesh.  For instance,
we often use the `"surface`" domain and the `"snow`" domain as aliases
for the same mesh, but making the distinction allows better names
throughout the code and user input and output
(e.g. `"snow-evaporation`" vs `"surface-evaporation`").  Note the use
of the `"-`" delimiter to separate domains from variables.

Occassionally it is useful to have a runtime-determined set of
domains.  For instance, if one would like to use an extruded 3D mesh,
and then solve equations only on vertical columns of that mesh, one
might **index** the columns according to how many surface faces there
are in the mesh, and create a domain for each face.  Sets of domains
like this are called **domain sets** and denoted by
`"DOMAIN_SET:INDEX`", e.g. `"column:3`", or `"surface_column:0`".
Indices may also be an enumerated list of strings,
e.g. `"watershed:upstream`" and `"watershed:downstream`" may split a
watershed into two subdomains.
