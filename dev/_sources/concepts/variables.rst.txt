Variable Naming Conventions
===========================

Given these concepts, we now can define Amanzi-ATS's full naming
convention.  A single piece of data is always named according to the
following convention:

    DOMAIN-ROOT_NAME\@TAG

As the word **variable** is loosely used throughout our documentation
as either the root variable name (e.g. `"evaporation`") or as the
domain-specific name (e.g. `"surface-evaporation`") or even as the
data name (e.g. `"surface-evaporation@NEXT`"), we leave that as
context specific and instead refer to the **root name** here.

Note that, by default, an empty domain name is used to refer to the
default, background mesh which defines the entirety of area to be used
in the simulation.  Typically, this is the subsurface, 3D domain, but
not always.  If no domain name appears, it can be assumed that the
quantity is defined on the default domain, e.g. `"saturation_liquid`"
is a valid name defined on the default domain.

Additionally, it can be useful to refer to derivatives.  Derivatives
are named with the character `"d`" to indicate either a total or
partial derivative, and a `"|`" to separate the numerator from the
denominator.  So, `"dsaturation_liquid@NEXT|dpressure@NEXT`" defines
the total derivative of liquid saturation with respect to pressure, on
the default domain, at the `"NEXT`" tag.

Sometimes, particularly in checkpoint files, we need to determine
between the same variable, on the same domain, but on different
**entities** (e.g. node, edge, face, or cell) of the discrete mesh.
In such cases, we append a `".`" and the entity name,
e.g. `"darcy_flux.face`" or `"surface-pressure.cell`".  Also,
potentially we must determine between different **degrees of
freedom**, e.g. the x-component of the velocity field, or the H+
aqueous specie's concentration.  These are additionally appended to
variable names with a `".`",
e.g. `"surface-total_component_concentration.cell.H+`" or
`"velocity.cell.2`" (the z-component of the velocity field, as we
always start numbering with 0).

This completes the list of special characters, which cannot be used in
any other name:

  * `"-`" (dash) is used between domain names and variable names, e.g. `"surface-ponded_depth`"
  * `":`" (colon) is used in domain sets, between the domain set name and the index, e.g. `"column:0`" or `"watershed:upstream`"
  * `"|`" (pipe) is used in derivative names, e.g. `"dwater_content|dpressure`"
  * `".`" (period) is used in a fully qualified name, such as for components of a vector, e.g. `"saturation.cell.0`" or `"free_ion_concentration.cell.NO3-`"

Blank spaces should be avoided, but names with spaces are still valid
and do not break the code (though they may break visualization and
other scripts, where we automatically replace them with
underscores).  Prefer to use `"_`" (underscore) instead, which may
appear in any part of a name.

In general, Amanzi-ATS names are relatively verbose. This is
intentional -- we believe that verbose names should be preferred to
jargon or abbreviated names, as it makes output more self-documenting.

